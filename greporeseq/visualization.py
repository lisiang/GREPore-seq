# -*- coding: utf-8 -*-
# @Time    : 2021/4/6
# @Author  : Siang_Li
# @File    : visualization.py
# @Software: PyCharm
"""Generating files that IGV needed"""
import os
import subprocess
import logging

logger = logging.getLogger('root')
logger.propagate = False

def visualizing(reference, demultis, outpath):
    """Generate .sorted.bam and .sorted.bam.bai files"""
    temoutfile = f'{os.path.splitext(reference)[0]}.mmi'
    mmifile = f'{os.path.join(outpath, os.path.split(temoutfile)[1])}'
    make_cmd = f'minimap2 -x map-ont -d "{mmifile}" "{reference}"'
    subprocess.run(make_cmd, shell=True)

    demultiname = f'{os.path.split(demultis)[1]}'
    purename = demultiname.strip(".fastq.gz")
    outsam = f'{os.path.join(outpath, purename)}.sam'
    outbam = f'{os.path.join(outpath, purename)}.bam'
    outsorted_bam = f'{os.path.join(outpath, purename)}.sorted.bam'

    align_cmd = f'minimap2 -ax map-ont "{mmifile}" "{demultis}" > "{outsam}"'
    sam2bam_cmd = f'samtools view -bS "{outsam}" > "{outbam}"'
    sortedbam_cmd = f'samtools sort -O bam -o "{outsorted_bam}" -T temp "{outbam}"'
    index_cmd = f'samtools index "{outsorted_bam}"'

    subprocess.run(align_cmd, shell=True)
    subprocess.run(sam2bam_cmd, shell=True)
    subprocess.run(sortedbam_cmd, shell=True)
    subprocess.run(index_cmd, shell=True)
    os.remove(outsam)
    os.remove(outbam)
    return outsorted_bam
