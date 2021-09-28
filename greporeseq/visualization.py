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

temp = os.path.join(os.getcwd(), '.temp')


def make_mmi_file(fa_list, outpath):
    """make mmi files that used in minimap2"""
    logger.info('Generating mmi file...')
    mmi_outfile = []
    with open(temp, 'w') as t:
        for filename in fa_list:
            temoutfile = f'{os.path.splitext(filename)[0]}.mmi'
            outfile = f'{os.path.join(outpath, os.path.split(temoutfile)[1])}'
            mmi_outfile.append(outfile)
            make_cmd = f'minimap2 -x map-ont -d "{outfile}" "{filename}"'
            subprocess.run(make_cmd, shell=True, stdout=t)
    logger.info('mmi file has been generated')
    return mmi_outfile


def visualizing(mmi_list, demulti_list, outpath):
    """Generate .sorted.bam and .sorted.bam.bai files"""
    sorted_bam_list = []
    with open(temp, 'w') as t:
        for mmifile in mmi_list:
            sitename = f'{os.path.split(mmifile)[1]}'
            sitename = f'{sitename.strip(".mmi")}'
            i = 1
            for demultifile in demulti_list:
                # if sitename in demultifile:
                logger.info('Generating files needed for IGV')
                demultiname = f'{os.path.split(demultifile)[1]}'
                purename = demultiname.strip(".fastq.gz")
                outsam = f'{os.path.join(outpath, purename)}-{sitename}.sam'
                outbam = f'{os.path.join(outpath, purename)}-{sitename}.bam'
                outsorted_bam = f'{os.path.join(outpath, purename)}-{sitename}.sorted.bam'
                sorted_bam_list.append(outsorted_bam)
                align_cmd = f'minimap2 -ax map-ont "{mmifile}" "{demultifile}" > "{outsam}"'
                sam2bam_cmd = f'samtools view -bS "{outsam}" > "{outbam}"'
                sortedbam_cmd = f'samtools sort -O bam -o "{outsorted_bam}" -T temp "{outbam}"'
                index_cmd = f'samtools index "{outsorted_bam}"'
                i += 1
                try:
                    subprocess.run(align_cmd, shell=True, stdout=t)
                except Exception as e:
                    logger.error(e)
                else:
                    subprocess.run(sam2bam_cmd, shell=True, stdout=t)
                    subprocess.run(sortedbam_cmd, shell=True, stdout=t)
                    subprocess.run(index_cmd, shell=True, stdout=t)
                    print('\n')
                    os.remove(outsam)
                    os.remove(outbam)
    return sorted_bam_list
