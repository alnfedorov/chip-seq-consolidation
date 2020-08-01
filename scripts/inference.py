import numpy as np
import torch

import utils
from data.pod import IntervalReads
from tqdm import tqdm


def predict(model: torch.nn.Module, crop: torch.Tensor):
    with torch.no_grad():
        cropenrich, croppeaks = model(crop)
        cropenrich, croppeaks = cropenrich[0, 0].cpu(), \
                                croppeaks[0, 1].exp().cpu()
    return cropenrich, croppeaks


FRAGMENT_SIZE = 250
EFFECTIVE_GENOME_SIZE = 2.7e9

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
MODEL_PATH = "/data/model-final.pth"
MODEL = torch.load(MODEL_PATH).to(DEVICE).eval()

for experiment in ("ENCSR652QNW", "ENCSR930HLX", "ENCSR956CTX"):
    for preset in ("1m", "5m", "10m", "25m", "original"):
# preset = "5m"
# experiment = "ENCSR361FWQ"
        print()
        print(f"EXPERIMENT: {experiment}, PRESET: {preset}")
        READS_DIR = f"/data/encode/H3K4me3/{experiment}/simulated/subsample/{preset}/filtered-bam/"
        if preset == "original":
            READS_DIR = READS_DIR.replace(f"/simulated/subsample/{preset}", "/original")
        import os
        TREATMENT = [os.path.join(READS_DIR, file) for file in os.listdir(READS_DIR) if "control" not in file and file.endswith(".bam")]
        CONTROL = [os.path.join(READS_DIR, file) for file in os.listdir(READS_DIR) if "control" in file and file.endswith(".bam")]

        # if preset == "original":
        #     TREATMENT = tuple(x.replace(f"/simulated/subsample/{preset}", "/original") for x in TREATMENT)
        #     CONTROL = tuple(x.replace(f"/simulated/subsample/{preset}", "/original") for x in CONTROL)

        SAVETO_PEAKS_FILE = f"/data/encode/viz/{experiment}-{preset}-peaks.bed"
        # SAVETO_PEAKS_SCORE_FILE = f"/data/encode/viz/{experiment}-{preset}-peak-scores.bdg"
        # SAVETO_SOURCE_ENRICHMENT_FILE = f"/data/encode/viz/{experiment}-{preset}-enrichment.bdg"
        # SAVETO_PREDICT_ENRICHMENT_FILE = f"/data/encode/viz/{experiment}-{preset}-pred-enrichment.bdg"

        BIN_SIZE = 25
        WINDOW_SIZE = 1024
        CROP_WINDOWS = 20
        OVERLAP_WINDOWS = 1

        MIN_PEAK_SIZE = 200

        treatment, chromosomes = utils.cached.separate_reads(TREATMENT)
        control, cchromosomes = utils.cached.separate_reads(CONTROL)
        assert all(chromosomes[k] == cchromosomes[k] for k in chromosomes) and \
               set(chromosomes.keys()) == set(cchromosomes.keys())

        SAVETO_PEAKS = open(SAVETO_PEAKS_FILE, 'w')
        # SAVETO_PEAKS_SCORE = open(SAVETO_PEAKS_SCORE_FILE, 'w')
        # SAVETO_SOURCE_ENRICHMENT = open(SAVETO_SOURCE_ENRICHMENT_FILE, 'w')
        # SAVETO_PREDICT_ENRICHMENT = open(SAVETO_PREDICT_ENRICHMENT_FILE, 'w')
        #
        # for file in (SAVETO_PEAKS_SCORE, SAVETO_PREDICT_ENRICHMENT, SAVETO_SOURCE_ENRICHMENT):
        #     file.write("track type=bedGraph\n")

        for chr in sorted(chromosomes):
            if chr.lower() == 'chrm':
                continue
            print(f"Start processing {chr}")
            chromlen = chromosomes[chr]

            treads = IntervalReads(np.load(treatment.forward[chr]), np.load(treatment.reverse[chr]))
            creads = IntervalReads(np.load(control.forward[chr]), np.load(control.reverse[chr]))
            if treads.numreads == 0 or creads.numreads == 0:
                continue

            enrich_ends, enrich_values = utils.enrichment.endtoend(chromlen, FRAGMENT_SIZE, EFFECTIVE_GENOME_SIZE,
                                                   treatment.numreads, control.numreads,
                                                   treads.forward, treads.reverse,
                                                   creads.forward, creads.reverse)
            enrichment = utils.enrichment.todense(enrich_ends, enrich_values, 0, chromlen)
            assert enrichment.size == chromlen

            # save enrichment
            # curstart, curend = 0, enrich_ends[0]
            # curval = round(enrich_values[0], 3)
            # for end, val in tqdm(zip(enrich_ends, enrich_values), total=enrich_ends.size):
            #     val = round(val, 3)
            #     if val != curval:
            #         SAVETO_SOURCE_ENRICHMENT.write(f"{chr}\t{curstart}\t{curend}\t{curval:.3f}\n")
            #         curstart, curend, curval = curend, end, val
            #     else:
            #         curend = end
            # SAVETO_SOURCE_ENRICHMENT.write(f"{chr}\t{curstart}\t{curend}\t{curval:.3f}")


            binned_enrichment = enrichment[:-(enrichment.size % BIN_SIZE)].reshape(-1, BIN_SIZE).mean(axis=-1)
            total_crops = len(binned_enrichment) // (CROP_WINDOWS * WINDOW_SIZE)
            assert total_crops > 0
            binned_enrichment = torch.tensor(binned_enrichment)

            end = CROP_WINDOWS * WINDOW_SIZE
            crop = binned_enrichment[:end].reshape(1, 1, -1).to(DEVICE)
            cropenrich, croppeaks = predict(MODEL, crop)
            predpeaks, predenrichment = [croppeaks], [cropenrich]

            overlap = OVERLAP_WINDOWS * WINDOW_SIZE
            for ind in range(1, total_crops):
                start, end = (CROP_WINDOWS * ind - OVERLAP_WINDOWS) * WINDOW_SIZE, \
                             CROP_WINDOWS * (ind + 1) * WINDOW_SIZE
                crop = binned_enrichment[start: end].reshape(1, 1, -1).to(DEVICE)

                cropenrich, croppeaks = predict(MODEL, crop)
                predpeaks[-1][-overlap:] = (predpeaks[-1][-overlap:] + croppeaks[:overlap]) / 2.0
                predenrichment[-1][-overlap:] = (predenrichment[-1][-overlap:] + cropenrich[:overlap]) / 2.0

                predpeaks.append(croppeaks[overlap:])
                predenrichment.append(cropenrich[overlap:])

            left = len(binned_enrichment) - total_crops * (CROP_WINDOWS * WINDOW_SIZE)
            if left != 0:
                padding = 256 - left % 256
                overlap = padding + OVERLAP_WINDOWS * WINDOW_SIZE
                assert overlap < CROP_WINDOWS * WINDOW_SIZE

                crop = binned_enrichment[-(left + overlap):].reshape(1, 1, -1).to(DEVICE)
                assert crop.shape[-1] % 256 == 0

                croppeaks, cropenrich = predict(MODEL, crop)
                predpeaks[-1][-overlap:] = (predpeaks[-1][-overlap:] + croppeaks[:overlap]) / 2.0
                predenrichment[-1][-overlap:] = (predenrichment[-1][-overlap:] + cropenrich[:overlap]) / 2.0

                predpeaks.append(croppeaks[overlap:])
                predenrichment.append(cropenrich[overlap:])

            assert sum(x.numel() for x in predpeaks) == sum(x.numel() for x in predenrichment) == binned_enrichment.numel()
            predpeaks, predenrichment = torch.cat(predpeaks).numpy(), torch.cat(predenrichment).numpy()

            # 1. Write peaks to the bed file
            # TODO: rewrite with numba(pure loop, no diff)
            tmp = np.hstack([[False], predpeaks > 0.5, [False]])  # padding
            diff = np.diff(tmp.astype(np.int8))
            starts, ends = np.where(diff == 1)[0], \
                           np.where(diff == -1)[0]

            for start, end in tqdm(zip(starts, ends), total=starts.size):
                SAVETO_PEAKS.write(
                    f"{chr}\t{start * BIN_SIZE}\t{end * BIN_SIZE}\t.\t0\t.\t{predpeaks[start:end].mean():.5f}\n"
                )

            # # 2. Write peak intensities and predicted enrichment
            # for file, values in zip([SAVETO_PEAKS_SCORE, SAVETO_PREDICT_ENRICHMENT],
            #                         [predpeaks, predenrichment]):
            #     curstart = 0
            #     curval = round(values[0], 3)
            #     for binind, val in tqdm(enumerate(values), total=values.size):
            #         val = round(val, 3)
            #         if val != curval:
            #             curstart, curend = curstart * BIN_SIZE, binind * BIN_SIZE
            #             file.write(f"{chr}\t{curstart}\t{curend}\t{curval:.3f}\n")
            #             curstart = binind
            #             curval = val
            #
            #     curstart, curend = curstart * BIN_SIZE, (values.size - 1) * BIN_SIZE + 1
            #     file.write(f"{chr}\t{curstart}\t{curend}\t{curval:.3f}")

        SAVETO_PEAKS.close()
        # SAVETO_PEAKS_SCORE.close()
        # SAVETO_SOURCE_ENRICHMENT.close()
        # SAVETO_PREDICT_ENRICHMENT.close()

#
# async def make_viewable():
#     import asyncio
#     from asyncio.subprocess import create_subprocess_shell
#     tool = "/data/encode/bedGraphToBigWig"
#     chrominfo = "/data/hg19/chromInfo.txt"
#
#     subprocesses = [
#         asyncio.create_task(create_subprocess_shell(f"{tool} {file} {chrominfo} {file}.bw"))
#         for file in [SAVETO_PREDICT_ENRICHMENT_FILE, SAVETO_SOURCE_ENRICHMENT_FILE, SAVETO_PEAKS_SCORE_FILE]
#     ]
#     await asyncio.gather(*subprocesses)
#
#     subprocesses = [
#         asyncio.create_task(create_subprocess_shell(f"rm {file}"))
#         for file in [SAVETO_PREDICT_ENRICHMENT_FILE, SAVETO_SOURCE_ENRICHMENT_FILE, SAVETO_PEAKS_SCORE_FILE]
#     ]
#     await asyncio.gather(*subprocesses)
#
# import asyncio
# asyncio.run(make_viewable())
#
#
# if preset == "original":
#     # apt install tabix
#     from data.parse import fromjson
#     consensus = fromjson("/data/encode/H3K4me3/ENCSR361FWQ/")[0].peaks
#     import os
#     import shutil
#     root = "/data/encode/viz"
#     shutil.move(consensus, os.path.join(root, 'consensus.bed'))
#
#     tozipindex = [x for x in os.listdir(root) if x.endswith(".bed")]
#     tozipindex = [os.path.join(root, x) for x in tozipindex] + [os.path.join(root, 'consensus.bed')]
#
#     for file in tozipindex:
#         saveto = f"{root}/{os.path.basename(file)}.gz"
#         if os.path.exists(saveto):
#             os.remove(saveto)
#             assert not os.path.exists(saveto), saveto
#         os.system(f'bgzip --stdout -@ 12 "{file}" > "{saveto}"')
#         os.system(f'tabix -p bed -f "{saveto}"')
#
#     for file in tozipindex:
#         if os.path.exists(file):
#             os.remove(file)
