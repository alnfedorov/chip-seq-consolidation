import numba
import numpy as np
from numba.core.types import int32, int64, float32


@numba.jit(cache=True, nopython=True)
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
    if ends_of_intervals[real_len-1] < rightmost_coord:
        ends_of_intervals[real_len] = rightmost_coord
        values[real_len] = 0
        real_len += 1

    assert ends_of_intervals[real_len-1] == rightmost_coord
    return ends_of_intervals[:real_len], values[:real_len]


def pileup(forward_reads_starts: np.ndarray, reverse_reads_starts: np.ndarray,
           five_shift: int32, three_shift: int32,
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
    :param five_shift: Optional 5' shift for all reads.
    :param three_shift: Optional 3' shift for all reads. Fragment size / 2, for example.
    :param leftmost_coord: inclusive, min position value, used to clip reads exceeding chromosome size
    :param rightmost_coord: exclusive, max position value, used to clip reads exceeding chromosome size
    :param scale: all pileup values will be multiplied by this scale
    :param baseline: baseline pileup value
    """
    total_reads = forward_reads_starts.size + reverse_reads_starts.size
    start_coords, end_coords = np.empty(total_reads, np.int32), np.empty(total_reads, np.int32)

    freads = forward_reads_starts.size
    # handle forward reads
    start_coords[:freads] = forward_reads_starts - five_shift
    end_coords[:freads] = forward_reads_starts + three_shift

    # handle reverse reads
    start_coords[freads:] = reverse_reads_starts - three_shift
    end_coords[freads:] = reverse_reads_starts + five_shift

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


@numba.jit(cache=True, nopython=True)
def per_interval_max(ends_of_intervals, values, scale: float32, baseline: float32):
    for p, v in zip(ends_of_intervals, values):
        assert p.size == v.size

    # 1. Make buffers for current intervals and results
    tracks = len(ends_of_intervals)
    track_nextind = np.zeros(tracks, dtype=np.int32)
    track_curend = np.empty(tracks, dtype=np.int32)
    track_curval = np.empty(tracks, dtype=np.float32)

    for track in range(tracks):
        track_curend[track] = ends_of_intervals[track][0]
        track_curval[track] = values[track][0]

    maxlength = int32(0)
    # sum(x.size for x in positions)
    for p in ends_of_intervals:
        maxlength += p.size
    res_ends = np.empty(maxlength, dtype=np.int32)
    res_values = np.empty(maxlength, dtype=np.float32)

    curval = max(baseline, np.max(track_curval))
    curend = min(track_curend)
    real_length = 0

    # ahead declaration of buffer to avoid memory allocation on each iteration
    to_drop = np.empty(tracks, dtype=np.int32)
    while tracks > 0:
        # 1. Next value is a max among present intervals and baseline value
        nextval = max(baseline, np.max(track_curval))

        # 2. End of the interval is a min among active intervals ends
        nextend = min(track_curend)

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
            if track_curend[track] == nextend:
                nextind = track_nextind[track] + 1
                # finished interval
                if nextind == ends_of_intervals[track].size:
                    to_drop[to_drop_total] = track
                    to_drop_total += 1
                    continue
                track_curend[track] = ends_of_intervals[track][nextind]
                track_curval[track] = values[track][nextind]
                track_nextind[track] = nextind

        if to_drop_total != 0:
            tracks -= to_drop_total
            indices = to_drop[:to_drop_total]
            track_curend = np.delete(track_curend, indices)
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


@numba.jit(cache=True, nopython=True)
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
