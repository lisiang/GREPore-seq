# -*- coding: utf-8 -*-
# @Time    : 2021/4/6
# @Author  : Siang_Li
# @File    : stats.py
# @Software: PyCharm
"""stating the input file's reads info"""
import os
import subprocess
import logging

logger = logging.getLogger('root')
logger.propagate = False

temp = os.path.join(os.getcwd(), '.temp')


def stats(demulted_path, output_path):
    logger.info('stating file reads info..')
    stat_output_file = []
    demulted_list = [os.path.join(demulted_path, file) for file in os.listdir(demulted_path)]
    for demulted_file in demulted_list:
        outfile = f'{os.path.join(output_path, "stat")}-result.txt'
        stat_output_file.append(outfile)
        statcmd = f'seqkit stat "{demulted_file}" -T >> "{outfile}"'
        with open(temp, 'w') as t:
            subprocess.call(statcmd, shell=True, stdout=t)

    for file in stat_output_file:
        con = []
        fr = open(file, 'r')
        for line in fr.readlines():
            con.append(line)
        context = con
        content = sorted(set(context), key=context.index)
        fw = open(file, 'w')
        fw.write(''.join(content))
        fw.close()


def coverage(visul_path, outpath):
    logger.info('calculating coverage info..')
    outlist = []
    mmilist = []
    sortedbamlist = []
    sitename_list = []
    for root, dircts, files in os.walk(visul_path):
        for file in files:
            if os.path.splitext(file)[1] == '.mmi':
                mmilist.append(os.path.join(visul_path, file))
                sitename_list.append(os.path.split(os.path.splitext(file)[0])[1])
            if os.path.splitext(file)[1] == '.bam':
                sortedbamlist.append(os.path.join(visul_path, file))

    with open(temp, 'w') as t:
        for sitename in sitename_list:
            for sorted_bam in sortedbamlist:
                if sitename in sorted_bam:
                    outname = os.path.join(outpath, f'{sitename}-coverage-result.txt')
                    outlist.append(outname)
                    count_cmd = f'samtools coverage "{sorted_bam}" >> "{outname}"'
                    subprocess.call(count_cmd, shell=True, stdout=t)

    for out in outlist:
        con = []
        fr = open(out, 'r')
        for line in fr.readlines():
            con.append(line)
        context = con
        content = sorted(set(context), key=context.index)
        fw = open(out, 'w')
        fw.write(''.join(content))
        fw.close()
