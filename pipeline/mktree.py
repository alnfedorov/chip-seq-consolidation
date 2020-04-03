import os
from .config import ORIGINAL_DIR, UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, UNIQUE_READS_DIR, PEAK_CALLING_DIR, \
    LOGS_DIR, CHIPS_DIR, CHIPS_MODELS_DIR, SIMULATED_DIR, SIMULATED_READS_DIRS, BAM_QC_DIR


def mktree(root: str):
    """
    Setup directories tree for the given experiment root
    Overall structure looks like this:
        ── original
        │   ├── duplicates-bam
        │   ├── logs
        │   ├── qc
        │   ├── peak-calling
        │   ├── unfiltered-bam
        │   └── unique-bam
        └── simulated
            └── chips
                ├── models
                ├── q0.25
                │   ├── duplicates-bam
                │   ├── logs
                │   ├── qc
                │   ├── peak-calling
                │   ├── unfiltered-bam
                │   └── unique-bam
                ├── q0.5
                .......
    :param root: path to the experiment root
    """
    assert os.path.isdir(root)

    def basic_tree(folder: str):
        for d in (LOGS_DIR, PEAK_CALLING_DIR, BAM_QC_DIR):
            d = os.path.join(folder, d)
            os.makedirs(d, exist_ok=True)
        for d in (UNFILTERED_READS_DIR, DUPLICATED_READS_DIR, UNIQUE_READS_DIR):
            d = os.path.join(folder, d)
            os.makedirs(d, exist_ok=True)

    original = os.path.join(root, ORIGINAL_DIR)
    basic_tree(original)

    chips = os.path.join(root, SIMULATED_DIR, CHIPS_DIR)
    for d in SIMULATED_READS_DIRS:
        d = os.path.join(chips, d)
        basic_tree(d)
    chips_params = os.path.join(chips, CHIPS_MODELS_DIR)
    os.makedirs(chips_params, exist_ok=True)
