# -*- coding: utf-8 -*-
# @Time         : 2021/8/6
# @Author       : Siang_Li
# @File         : condition_grep.py
# @Software     : PyCharm
# @E-mail       : lisiang.darcy@gmail.com
# @Description  : Grab read containing the specified number of sequences
import os
import shutil
import gzip
import itertools
import re
from collections import Counter


def set_dir(wd, direct_name, replace=True):
    """函数用于创建文件夹"""
    direct_abspath = os.path.join(wd, direct_name)
    if direct_name not in os.listdir(wd):
        os.mkdir(direct_abspath)
    elif replace:
        shutil.rmtree(direct_abspath)
        os.mkdir(direct_abspath)
    return direct_abspath


def read_fastq(fastq_name):
    fastqlist = []
    with gzip.open(fastq_name, 'rb') as fastq:
        while fastq:
            fastqinfo = {}
            fastqinfo['id'] = fastq.readline().decode().strip('\n')
            if not fastqinfo['id'] and fastq:
                break
            fastqinfo['seq'] = fastq.readline().decode().strip('\n')
            fastqinfo['direction'] = fastq.readline().decode().strip('\n')
            fastqinfo['quality'] = fastq.readline().decode().strip('\n')
            fastqlist.append(fastqinfo)
    return fastqlist


def grep(condition_grepseqs, fastqs, out_path, number=2):
    out_list = []
    for fastq in fastqs:
        fastq_name = os.path.split(fastq)[1]
        fastq_list = (read for read in read_fastq(fastq))
        for grepseq in condition_grepseqs:
            grepseq_name = os.path.split(grepseq)[1]
            readsid, pass_readsids = [], []
            output_name = os.path.join(out_path,
                                       f'{fastq_name.strip("fastq.gz")}-{grepseq_name.strip(".txt")}.fastq.gz')
            with open(grepseq, 'r') as f:
                sequence = [seq.strip() for seq in f.readlines()]
                for i in range(len(sequence)):
                    fastq_list, copy_fastqs = itertools.tee(fastq_list)
                    for read in copy_fastqs:
                        if re.search(sequence[i], read['seq']):
                            readsid.append(read['id'].strip().split()[0][1:])

            dup_dic = Counter(readsid)
            for k, v in dup_dic.items():
                if v >= number:
                    pass_readsids.append(k)

            with gzip.open(output_name, 'wb') as outf:
                fastq_list, copy_fastqs2 = itertools.tee(fastq_list)
                n = 0
                for read2 in copy_fastqs2:
                    if read2['id'].strip().split()[0][1:] in pass_readsids:
                        n += 1
                        if n == 1:
                            outf.write(
                                f"{read2['id']}\n{read2['seq']}\n{read2['direction']}\n{read2['quality']}".encode())
                            continue
                        text = f"\n{read2['id']}\n{read2['seq']}\n{read2['direction']}\n{read2['quality']}"
                        outf.write(text.encode())
            out_list.append(output_name)
    return out_list