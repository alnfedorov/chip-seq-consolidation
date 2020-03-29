
# files must be sorted by the read names !!!!!!
# samtools sort -n sample.bam sample.sorted_by_name


def run():
#  AT LEAST 2 replicas
#  might be single control for all of them
"""
PePr -c /data/encode/H3K4me3/ENCSR661MUS/data/duplicates-bam/duplicates-breplica-ENCFF755PXB.bam,/data/encode/H3K4me3/ENCSR661MUS/data/duplicates-bam/duplicates-control-ENCFF004LNU.bam \
     -i /data/encode/H3K4me3/ENCSR661MUS/data/duplicates-bam/duplicates-control-ENCFF004LNU.bam \
     -f bam \
     --threshold=0.00001 \
     --peaktype=sharp \
     --num-processors=12
"""
    # PePr postprocessing - skip for now
    pass