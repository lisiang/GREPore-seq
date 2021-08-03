# -*- coding: utf-8 -*-
# @Time    : 2021/4/5
# @Author  : Siang_Li
# @File    : demultiplex.py
# @Software: PyCharm
"""demultiplex functions"""
import os
import subprocess
import logging

logger = logging.getLogger('root')
logger.propagate = False

temp = os.path.join(os.getcwd(), '.temp')


def isolate_PCR_products(np_fr_data, grep_seqL, grep_seqR, outpath):
    """Isolate the PCR products"""
    logger.info('Isolating PCR products...')
    batch_number = os.path.split(np_fr_data)[1].strip('-chop-FR.fastq.gz')
    PCR_products = []
    i = 0
    while i < len(grep_seqL):
        site_name = [name[:name.index('GrepL')] for name in [os.path.split(n)[1] for n in grep_seqL]]
        out_name = os.path.join(outpath, f'{batch_number}-{site_name[i]}.fastq.gz')
        isolate_cmd = 'seqkit grep -s -f "{0}" "{1}" | seqkit grep -s -f "{2}" -o "{3}"'.format(grep_seqL[i],
                                                                                                np_fr_data,
                                                                                                grep_seqR[i],
                                                                                                out_name)
        with open(temp, 'w') as t:
            subprocess.call(isolate_cmd, shell=True, stdout=t)
        PCR_products.append(out_name)
        i += 1
    return PCR_products


def demultiplex_PCR_products(PCR_products, BCPSlist, outpath):
    """demultiplex the isolated PCR products"""
    logger.info('According BC-primer-seq demultiplex PCR products...')
    demulted_products = []
    site_demulted_dict = {}
    for product in PCR_products:
        batch_site_info = os.path.split(product)[1].strip('.fastq.gz')
        sitename = batch_site_info[batch_site_info.rfind('-')+1:]
        i = 0
        while i < len(BCPSlist):
            if sitename in BCPSlist[i]:
                output_name = f'{os.path.join(outpath, os.path.splitext(os.path.split(BCPSlist[i])[1])[0])}.fq.gz'
                site_demulted_dict[output_name] = sitename
                demul_cmd = f'seqkit grep -s -R 1:20 -i -r -f "{BCPSlist[i]}" "{product}" -o "{output_name}"'
                demulted_products.append(output_name)
                with open(temp, 'w') as t:
                    subprocess.call(demul_cmd, shell=True, stdout=t)
            i += 1
    return demulted_products, site_demulted_dict
