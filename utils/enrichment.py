import numba
import numpy as np
from numba.core.types import int32, float32


@numba.jit(cache=True, nopython=True, nogil=True)
def _dopileup(start_coords, end_coords, rightmost_coord: int32, scale: float32, baseline: float32):
    # build pileup
    # positions -> start positions of the constant region
    # values -> constant region value
    real_len = 0
    # reserve maxlen aot to save time
    ends_of_intervals = np.empty(2 * start_coords.size + 1, dtype=np.int32)
    values = np.empty(2 * start_coords.size + 1, dtype=np.float32)

    start_idx, end_idx = 0, 0

    pileup = int32(0)
    prev_coord = min(start_coords[start_idx], end_coords[end_idx])

    if prev_coord != 0:
        ends_of_intervals[real_len] = prev_coord
        values[real_len] = baseline
        real_len += 1

    total_read = start_coords.size
    while start_idx < total_read and end_idx < total_read:
        if start_coords[start_idx] < end_coords[end_idx]:
            coord = start_coords[start_idx]
            if coord != prev_coord:
                ends_of_intervals[real_len] = coord
                values[real_len] = max(pileup * scale, baseline)
                real_len += 1
                prev_coord = coord
            pileup += 1
            start_idx += 1
        elif start_coords[start_idx] > end_coords[end_idx]:
            coord = end_coords[end_idx]
            if coord != prev_coord:
                ends_of_intervals[real_len] = coord
                values[real_len] = max(pileup * scale, baseline)
                real_len += 1
                prev_coord = coord
            pileup -= 1
            end_idx += 1
        else:
            start_idx += 1
            end_idx += 1

    if end_idx < total_read:
        for i in range(end_idx, total_read):
            coord = end_coords[i]
            if coord != prev_coord:
                ends_of_intervals[real_len] = coord
                values[real_len] = max(pileup * scale, baseline)
                prev_coord = coord
                real_len += 1
            pileup -= 1
    assert pileup == 0

    # force intervals to occupy all space
    if ends_of_intervals[real_len - 1] < rightmost_coord:
        ends_of_intervals[real_len] = rightmost_coord
        values[real_len] = 0
        real_len += 1

    assert ends_of_intervals[real_len - 1] == rightmost_coord
    return ends_of_intervals[:real_len], values[:real_len]


def pileup(forward_reads_starts: np.ndarray, reverse_reads_starts: np.ndarray,
           five_extension: int32, three_extension: int32,
           leftmost_coord: int32, rightmost_coord: int32, scale: float32, baseline: float32):
    """
    Compute interval-wise pileup values, i.e. for each interval it is a number of DNA fragments covering the interval
    Each interval pileup is max(baseline, pileup * scale).
    Returned intervals are half-open[) and represented by two arrays.
    1. ends_of_intervals - end positions of the intervals. Start positions of the intervals sequence is 0
    2. values - pileup values
    So, i-th interval has coordinates [positions[i-1], positions[i]) and value values[i]. First interval starts from 0.
    :param forward_reads_starts: start positions of the forward reads
    :param reverse_reads_starts: start positions of the reverse reads
    :param five_extension: Optional 5' extension for all reads.
    :param three_extension: Optional 3' extension for all reads. Fragment size / 2, for example.
    :param leftmost_coord: inclusive, min position value, used to clip reads exceeding chromosome size
    :param rightmost_coord: exclusive, max position value, used to clip reads exceeding chromosome size
    :param scale: all pileup values will be multiplied by this scale
    :param baseline: baseline pileup value
    """
    total_reads = forward_reads_starts.size + reverse_reads_starts.size
    start_coords, end_coords = np.empty(total_reads, np.int32), np.empty(total_reads, np.int32)

    freads = forward_reads_starts.size
    # handle forward reads
    start_coords[:freads] = forward_reads_starts - five_extension
    end_coords[:freads] = forward_reads_starts + three_extension

    # handle reverse reads
    start_coords[freads:] = reverse_reads_starts - three_extension
    end_coords[freads:] = reverse_reads_starts + five_extension

    # clip coordinates in-place
    start_coords = np.clip(start_coords, leftmost_coord, rightmost_coord - 1, out=start_coords)
    end_coords = np.clip(end_coords, leftmost_coord, rightmost_coord, out=end_coords)

    # sort for simplicity
    start_coords = np.sort(start_coords)
    end_coords = np.sort(end_coords)

    assert start_coords.shape == end_coords.shape
    # run compiled function to actually compute the coverage of dna by fragments
    ends_of_intervals, values = _dopileup(start_coords, end_coords, rightmost_coord, scale, baseline)
    assert ends_of_intervals[-1] == rightmost_coord
    return ends_of_intervals, values


@numba.jit("float32[::1](int32[::1], float32[::1], int32, int32)", cache=True, nopython=True, nogil=True)
def todense(ends_of_intervas, values, dense_start, dense_end):
    assert ends_of_intervas.size == values.size
    if ends_of_intervas.size == 0:
        return np.zeros(dense_end - dense_start, dtype=np.float32)

    result = np.empty(dense_end - dense_start, dtype=np.float32)
    begin = dense_start
    for end, value in zip(ends_of_intervas, values):
        if end <= begin:
            continue
        s, e = begin - dense_start, end - dense_start
        result[s: e] = value
        begin = end
        if end >= dense_end:
            break
    assert end >= dense_end
    return result


@numba.jit(cache=True, nopython=True, nogil=True)
def per_interval_max(ends_of_intervals, values, scale: float32, baseline: float32):
    for p, v in zip(ends_of_intervals, values):
        assert p.size == v.size

    # 1. Make buffers for current intervals and results
    tracks = len(ends_of_intervals)
    track_nextind = np.zeros(tracks, dtype=np.int32)
    track_curent = np.empty(tracks, dtype=np.int32)
    track_curval = np.empty(tracks, dtype=np.float32)

    for track in range(tracks):
        track_curent[track] = ends_of_intervals[track][0]
        track_curval[track] = values[track][0]

    maxlength = int32(0)
    # sum(x.size for x in positions)
    for p in ends_of_intervals:
        maxlength += p.size
    res_ends = np.empty(maxlength, dtype=np.int32)
    res_values = np.empty(maxlength, dtype=np.float32)

    curval = max(baseline, np.max(track_curval))
    nextval = curval

    curend = min(track_curent)
    real_length = 0

    # ahead declaration of buffer to avoid memory allocation on each iteration
    to_drop = np.empty(tracks, dtype=np.int32)
    while tracks > 0:
        # 1. Next value is a max among present intervals and baseline value
        nextval = max(baseline, np.max(track_curval))

        # 2. End of the interval is a min among active intervals ends
        nextend = min(track_curent)

        # 3. Save interval if new value is encountered
        if nextval != curval:
            res_values[real_length] = curval * scale
            res_ends[real_length] = curend
            real_length += 1
            curval = nextval
        curend = nextend

        # 4. Push intervals if needed and drop finished intervals
        to_drop_total = 0
        for track in range(tracks):
            if track_curent[track] == nextend:
                nextind = track_nextind[track] + 1
                # finished interval
                if nextind == ends_of_intervals[track].size:
                    to_drop[to_drop_total] = track
                    to_drop_total += 1
                    continue
                track_curent[track] = ends_of_intervals[track][nextind]
                track_curval[track] = values[track][nextind]
                track_nextind[track] = nextind

                # update max value for the current interval
                # if values[track][nextind - 1] == nextval:
                #     nextval = max(baseline, np.max(track_curval))
                # else:
                #     nextval = max(nextval, track_curval[track])

        if to_drop_total != 0:
            tracks -= to_drop_total
            indices = to_drop[:to_drop_total]
            track_curent = np.delete(track_curent, indices)
            track_curval = np.delete(track_curval, indices)
            track_nextind = np.delete(track_nextind, indices)

            for track in to_drop[:to_drop_total][::-1]:
                assert track < len(ends_of_intervals) and track < len(values) and \
                       len(ends_of_intervals) == len(values) and \
                       len(ends_of_intervals) > 0
                ends_of_intervals.pop(track)
                values.pop(track)

    res_values[real_length] = curval * scale
    res_ends[real_length] = curend
    real_length += 1

    assert res_ends.size >= real_length
    return res_ends[:real_length], res_values[:real_length]


@numba.jit(cache=True, nopython=True, nogil=True)
def enrichment(treatment_endcoords, treatment_values, control_endcoords, control_values, pseudocount: float32):
    # assert treatment_endpos.size == treatment_val.size and control_endpos.size == control_val.size
    # assert treatment_endpos.dtype == control_endpos.dtype == np.int32 and \
    #        treatment_val.dtype == control_val.dtype == np.float32
    assert treatment_endcoords[-1] == control_endcoords[-1]

    maxlen = treatment_values.size + control_values.size
    res_ends = np.empty(maxlen, dtype=np.int32)
    res_values = np.empty(maxlen, dtype=np.float32)

    treatment_ind, control_ind = 0, 0
    tend, cend = treatment_endcoords[treatment_ind], control_endcoords[control_ind]
    tvalue, cvalue = treatment_values[treatment_ind], control_values[control_ind]

    curvalue = (tvalue + pseudocount) / (cvalue + pseudocount)
    curend = min(tend, cend)
    real_len = 0

    while True:
        # 1. Make next interval endpoint and value
        nextval = (tvalue + pseudocount) / (cvalue + pseudocount)
        nextend = min(tend, cend)

        # 2. Save if needed
        if nextval != curvalue:
            res_ends[real_len] = curend
            res_values[real_len] = curvalue
            real_len += 1
            curvalue = nextval
        curend = nextend

        # 3. Advance if needed
        if nextend == tend:
            treatment_ind += 1
            if treatment_ind < treatment_endcoords.size:
                tend = treatment_endcoords[treatment_ind]
                tvalue = treatment_values[treatment_ind]
        if nextend == cend:
            control_ind += 1
            if control_ind < control_endcoords.size:
                cend = control_endcoords[control_ind]
                cvalue = control_values[control_ind]

        if treatment_ind == treatment_endcoords.size or control_ind == control_endcoords.size:
            break

    res_ends[real_len] = curend
    res_values[real_len] = curvalue
    real_len += 1

    assert treatment_ind == treatment_values.size and control_ind == control_values.size
    assert res_ends[real_len - 1] == treatment_endcoords[-1] == control_endcoords[-1]
    return res_ends[:real_len], res_values[:real_len]


# @numba.jit(cache=True, nopython=True, nogil=True)
def endtoend(chromlen: int32, fragment_size: int32, effective_genome_size: float32,
             alltreatment_reads: int32, allcontrol_reads: int,
             treatment_forward: np.ndarray, treatment_reverse: np.ndarray,
             control_forward: np.ndarray, control_reverse: np.ndarray):
    """alltreatment_reads and allcontrol_reads ARE FROM THE WHOLE EXPERIMENT, NOT ONLY THIS CHR"""
    if treatment_forward.size == 0 and treatment_reverse.size == 0:
        return np.empty(0, dtype=np.int32), np.empty(0, dtype=np.float32)

    # 1. treatment pileup
    treatment_ends, treatment_pileup = pileup(treatment_forward, treatment_reverse,
                                              five_extension=0, three_extension=fragment_size,
                                              leftmost_coord=0, rightmost_coord=chromlen,
                                              scale=1.0, baseline=0.0)
    # 2. control lambda
    half_fragment = fragment_size // 2
    cends, cpileup = pileup(control_forward, control_reverse,
                            five_extension=half_fragment, three_extension=half_fragment,
                            leftmost_coord=0, rightmost_coord=chromlen, scale=1.0, baseline=0.0)
    cends1k, clambda1k = pileup(control_forward, control_reverse,
                                five_extension=500, three_extension=500,
                                leftmost_coord=0, rightmost_coord=chromlen,
                                scale=fragment_size / 1000, baseline=0.0)
    cends10k, clambda10k = pileup(control_forward, control_reverse,
                                  five_extension=5000, three_extension=5000,
                                  leftmost_coord=0, rightmost_coord=chromlen,
                                  scale=fragment_size / 10000, baseline=0.0)

    background_lambda = allcontrol_reads * fragment_size / effective_genome_size
    treatment_control_ratio = alltreatment_reads / allcontrol_reads

    cends, clambda = per_interval_max([cends, cends1k, cends10k], [cpileup, clambda1k, clambda10k],
                                      scale=treatment_control_ratio, baseline=background_lambda)

    # 3. enrichment
    enrich_ends, enrich_values = enrichment(treatment_ends, treatment_pileup, cends, clambda,
                                            pseudocount=0.0)
    return enrich_ends, enrich_values
