# -*- coding: utf-8 -*-
# @Time    : 2021/12/13 11:02
# @Author  : SiangLee
# @File    : read_fastq.py
# Description: Read fastq formate file
import gzip
import time

import logging

logger = logging.getLogger('root')
logger.propagate = False


def read_fastq(fastq):
    logger.info(f'Reading {fastq}')
    if fastq.endswith('.gz'):
        fq = gzip.open(fastq, 'rb')
    else:
        fq = open(fastq, 'r')
    with fq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            elif l1.decode('utf-8') == '\n':
                continue
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            read = l1 + l2 + l3 + l4
            if isinstance(read, bytes):
                read = str(read, encoding='utf-8')
            read = tuple(filter(None, read.split('\n')))
            yield read
