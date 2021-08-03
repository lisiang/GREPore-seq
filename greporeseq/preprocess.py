# -*- coding: utf-8 -*-
# @Time    : 2021/4/5
# @Author  : Siang_Li
# @File    : preprocess.py
# @Software: PyCharm
import logging
import os
import subprocess
import gzip
import time

import xlrd

logger = logging.getLogger('root')
logger.propagate = False

temp = os.path.join(os.getcwd(), '.temp')


def rm_dup(listA):
    return list(sorted(set(listA), key=listA.index))


def rm_none(listA):
    return list(filter(None, listA))


def read_excel(excel_name, cols, upper=False, dup_rm=True, none_rm=True):
    """read excel contents"""
    try:
        wb = xlrd.open_workbook(filename=excel_name)
        sheet = wb.sheet_by_index(0)
    except FileNotFoundError:
        logger.error(
            'No such file or directory: {0}. Please ensure your provided correct excel file.'.format(
                excel_name))
    else:
        contents = []
        i = 1
        while i < sheet.nrows:
            values = sheet.cell(i, cols).value
            if upper:
                contents.append(str(values).replace(' ', '').replace('\n', '').replace('\t', '').upper())
            else:
                contents.append(str(values).replace(' ', '').replace('\n', '').replace('\t', ''))
            i += 1
        if dup_rm:
            contents = rm_dup(contents)
        if none_rm:
            contents = rm_none(contents)
        return contents


def read_fastq(fastq):
    """read fastq file content"""
    if fastq.endswith('.gz'):
        fq = gzip.open(fastq, 'rb')
    else:
        fq = open(fastq, 'r')
    with fq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            read = l1 + l2 + l3 + l4
            if isinstance(read, bytes):
                read = str(read, encoding='utf-8')
            read = tuple(filter(None, read.split('\n')))
            yield read


def seq_rev_com(read):
    """reverse complement"""
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    rev_com_seq = read[1].translate(trantab)[::-1]
    line = f"\n{read[0]}\n{rev_com_seq}\n{read[2]}\n{read[3]}\n{read[0]}\n{read[1]}\n{read[2]}\n{read[3]}"
    return line


def rev_com(fastq):
    start = time.time()
    total_count = 1
    output = f"{fastq.replace('.fastq.gz', '-FR.fastq.gz')}"
    with gzip.open(output, 'wb')as f:
        for read in read_fastq(fastq):
            total_count += 1
            f.write(seq_rev_com(read).encode())
            if total_count % 100000 == 0:
                logger.info("Processed %d reads in %.1f minutes.", total_count, (time.time() - start) / 60)
    return output


def consolidate_np_data(np_data_name):
    """consolidation of ONT sequencing data"""
    if np_data_name.endswith('.fq.gz'):
        np_data_name = np_data_name.replace('.fq.gz', '.fastq.gz')

    if 'chop-fr' in np_data_name.lower():
        logger.info('Consolidation completed')
        return np_data_name

    elif 'chop' not in np_data_name.lower():
        logger.info('Chopping adapter...')
        choped_np_data = f"{np_data_name.replace('.fastq.gz', '-chop.fastq.gz')}"
        chop_cmd = f'porechop --adapter_threshold 85 --extra_end_trim 0 -i {np_data_name} -o "{choped_np_data}"'
        with open(temp, 'w') as t:
            subprocess.call(chop_cmd, shell=True, stdout=t)
        if '-fr' not in np_data_name.lower():
            logger.info('Merging forward and reverse...')
            choped_fr = rev_com(choped_np_data)
            return choped_fr
        return choped_np_data

    elif '-fr' not in np_data_name.lower():
        logger.info('Merging forward and reverse...')
        choped_fr_gz = rev_com(np_data_name)
        return choped_fr_gz


def make_fa(site_names, site_seqs, out_path):
    """generate the fa file"""
    logger.info('Generating AMP fa file...')
    txtoutfilelist = []
    fa_out_list = []
    i = 0
    while i < len(site_names):
        txt_out_file = os.path.join(out_path, f"{site_names[i]}.txt")
        txtoutfilelist.append(txt_out_file)
        file = open(txt_out_file, 'w')
        file.write('>' + site_names[i] + '\n' + site_seqs[i])
        file.close()

        path_site_name = txt_out_file.strip('.txt')
        txt2fasta_cmd = 'seqkit fx2tab "{0}" | seqkit tab2fx | seqkit seq -u | seqkit replace -p " |-" -s > "{1}.fasta"'.format(
            txt_out_file, path_site_name)
        with open(temp, 'w') as t:
            subprocess.call(txt2fasta_cmd, shell=True, stdout=t)
        fa_out_list.append(f'{path_site_name}.fasta')
        os.remove(txt_out_file)
        i += 1
    return fa_out_list


def make_bcps(bc_primer_names, bc_primer_seqs, bc_primer_lenths, output_path):
    """generate the BC-primer-seq file"""
    logger.info('Generating BC-primer-seq file...')
    BC_primer_seqs = []
    fa_list = make_fa(bc_primer_names, bc_primer_seqs, output_path)
    i = 0
    while i < len(bc_primer_names):
        with open(temp, 'w') as t:
            fa_txt_temp = os.path.join(output_path, f'{bc_primer_names[i]}.temp')
            makeBCP_cmd = f'seqkit subseq -r 1:{int(float(bc_primer_lenths[i]) + 5)} "{fa_list[i]}" | seqkit sliding -W 9 -s 1 | seqkit fx2tab >> "{fa_txt_temp}"'
            subprocess.call(makeBCP_cmd, shell=True, stdout=t)

            BCPout = os.path.join(output_path, f'{bc_primer_names[i]}.txt')
            with open(fa_txt_temp, 'r') as fatxttemp:
                f = open(BCPout, 'w')
                for line in fatxttemp.readlines():
                    line = line.split('\t')
                    f.write(line[1] + '\n')
                f.close()
            os.remove(fa_txt_temp)
            os.remove(fa_list[i])
            os.remove(fa_list[i] + '.seqkit.fai')
            BC_primer_seqs.append(BCPout)
        i += 1
    return BC_primer_seqs


def make_grepseq(amp_fa_list):
    """Generate Grepseq file by AMP fa file"""
    logger.info('Generating grepseq file...')
    GPseqlist = []
    for fapathname in amp_fa_list:
        sitename = fapathname[:fapathname.index('.fasta')]
        tmpL, tmpR = sitename + 'GrepL-15_5F.tmp', sitename + 'GrepR-15_5F.tmp'
        outL, outR = sitename + 'GrepL-15_5F.txt', sitename + 'GrepR-15_5F.txt'

        maketmpL = f'seqkit subseq -r 20:90 "{fapathname}" | seqkit sliding -W 15 -s 5 | seqkit fx2tab >> "{tmpL}"'
        maketmpR = f'seqkit subseq -r -90:-20 "{fapathname}" | seqkit sliding -W 15 -s 5 | seqkit fx2tab >> "{tmpR}"'

        with open(temp, 'w') as t:
            subprocess.call(maketmpL, shell=True, stdout=t)
            subprocess.call(maketmpR, shell=True, stdout=t)

        with open(tmpL, 'r', encoding='utf-8') as tempfL:
            with open(outL, 'w') as fL:
                for line in tempfL.readlines():
                    line = line.split('\t')
                    fL.write(line[1] + '\n')
                fL.close()
            GPseqlist.append(outL)

        with open(tmpR, 'r', encoding='utf-8') as tempfR:
            with open(outR, 'w') as fR:
                for line in tempfR.readlines():
                    line = line.split('\t')
                    fR.write(line[1] + '\n')
                fR.close()
            GPseqlist.append(outR)

        os.remove(fapathname + '.seqkit.fai')
        os.remove(tmpL)
        os.remove(tmpR)
    return GPseqlist
