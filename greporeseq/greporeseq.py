# -*- coding: utf-8 -*-
# @Time         : 2021/7/27
# @Author       : Siang_Li
# @File         : greporeseq_ori.py
# @Software     : PyCharm
# @E-mail       : lisiang.darcy@gmail.com
# @Description  : a pipeline for analysing ONT sequencing data

import os
import argparse
import re
import textwrap
import time
import subprocess
from tkinter.tix import Tree
import random

import log
import visualization as vis
import read_input as ri
from read_fastq import read_fastq

logger = log.createCustomLogger('root')


def rm_dup_none(l):
    return list(filter(None, sorted(set(l), key=l.index)))


def filter_none(l):
    return list(filter(None, l))


class GREPoreSeq:
    def __init__(self):
        pass

    def set_dir(self, dir_name):
        """Setting directory"""
        if dir_name not in os.listdir():
            os.mkdir(dir_name)
            return dir_name
        else:
            return dir_name

    def seq_rev_com(self, read):
        """reverse complement"""
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
        rev_com_seq = read[1].translate(trantab)[::-1]
        rev_com_id = f'{read[0].strip().split()[0]}_R'
        linel = [rev_com_id, rev_com_seq, read[2], read[3][::-1]]
        return linel

    def rev_com(self, seq):
        if seq:
            trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
            return seq.translate(trantab)[::-1]
        else:
            return

    def input_info(self, ref_yaml=None, demulti_yaml=None):
        """Extraction of excel information"""
        if ref_yaml:
            self.ref_info = ri.r_input(ref_yaml)
        if demulti_yaml:
            self.demul_info = ri.r_input(demulti_yaml)

    def mk_reference(self):
        """Check and create the output folder"""
        logger.info('Making FASTA file')
        self.ref_dir = self.set_dir('Reference')

        self.id_references = {}
        for ref_id, reference_info in self.ref_info.items():
            reference_id = self.reform_id(ref_id)
            seq = self.reform_seq(reference_info['sequence'])
            line = f'>{reference_id}\n'
            line += textwrap.fill(seq, 60)
            fa_out = os.path.join(self.ref_dir, f'{reference_id}.fasta')
            self.id_references[reference_id] = fa_out
            with open(fa_out, 'w') as f:
                f.write(line)

    def mk_grepseq(self, ref_seq, walk, step):
        if ref_seq:
            s, length = 0, len(ref_seq)
            greps = [ref_seq[s + step * i:walk + step * i] for i in range(round(length / step))]
            grepws = [seq for seq in greps if len(seq) == walk]
            return grepws
        else:
            return

    def osearch(slef, seqs_used, seq_searched, match):
        m = 0
        if seqs_used:
            for seq in seqs_used:
                if re.search(seq, seq_searched):
                    m += 1
                    if m >= match:
                        return m
        elif not seqs_used:
            return True

    def reform_seq(self, seq):
        if seq:
            return seq.replace(' ', '').replace('\n', '').replace('\t', '').upper()
        else:
            return

    def reform_num(self, num):
        if num:
            return int(num)
        else:
            return

    def reform_id(self, name_id):
        if name_id:
            return name_id.replace(' ', '_').replace('.', '_').replace('\t', '_')
        else:
            return

    def disassemble(self):
        logger.info('Disassembling input info')
        self.fastq_ids = []
        self.id_ref, self.left_seqs, self.right_seqs, self.BCprimer_Fs, self.BClen_Fs, self.id_BCprimer_Rs, self.id_uniseq, self.id_BClen_R = {}, {}, {}, {}, {}, {}, {}, {}
        for demuti_id, demuti_info in self.demul_info.items():
            demuti_id = self.reform_id(demuti_id)
            self.fastq_ids.append(demuti_id)
            self.left_seqs[f'{demuti_id}'] = self.reform_seq(demuti_info['left_seqs'])
            self.right_seqs[f'{demuti_id}'] = self.reform_seq(demuti_info['right_seqs'])
            self.BCprimer_Fs[f'{demuti_id}'] = self.reform_seq(demuti_info['BCprimer_F'])
            self.BClen_Fs[f'{demuti_id}'] = self.reform_num(demuti_info['BClen_F'])
            self.id_ref[f'{demuti_id}'] = self.reform_id(demuti_info['reference_id'])
            self.id_BCprimer_Rs[f'{demuti_id}'] = self.rev_com(self.reform_seq(demuti_info['BCprimer_R']))
            self.id_BClen_R[f'{demuti_id}'] = self.reform_num(demuti_info['BClen_R'])
            self.id_uniseq[f'{demuti_id}'] = self.reform_seq(demuti_info['unique_sequence'])

    def prepare_grepseq(self):
        logger.info('Preparing Grepseqs')
        self.id_reads_dic, self.id_grepseq_dic = {}, {}
        for f_id in self.fastq_ids:
            self.l_seq = self.left_seqs[f_id]
            self.r_seq = self.right_seqs[f_id]
            BCprimer_F = self.BCprimer_Fs[f_id]
            BCprimerR = self.id_BCprimer_Rs[f_id]
            BClen_F = self.BClen_Fs[f_id]
            BClen_R = self.id_BClen_R[f_id]
            uniseq = self.id_uniseq[f_id]

            if not BClen_F:
                BCrange_F = -1
            else:
                BCrange_F = BClen_F * 2 - 8

            if not BClen_R:
                BCrange_R = -1
            else:
                BCrange_R = BClen_R * 2 - 8

            self.id_reads_dic[f_id] = []
            if self.l_seq:
                self.id_grepseq_dic[f'{f_id}leftseq'] = self.mk_grepseq(ref_seq=self.l_seq, walk=17, step=20)
            else:
                self.id_grepseq_dic[f'{f_id}leftseq'] = None
            if self.r_seq:
                self.id_grepseq_dic[f'{f_id}rightseq'] = self.mk_grepseq(ref_seq=self.r_seq, walk=17, step=20)
            else:
                self.id_grepseq_dic[f'{f_id}rightseq'] = None
            if BCprimer_F:
                self.id_grepseq_dic[f'{f_id}BCprimer_F'] = self.mk_grepseq(ref_seq=BCprimer_F[:BCrange_F], walk=11,
                                                                           step=1)
            else:
                self.id_grepseq_dic[f'{f_id}BCprimer_F'] = None
            if BCprimerR:
                self.id_grepseq_dic[f'{f_id}BCprimer_R'] = self.mk_grepseq(ref_seq=BCprimerR[:BCrange_R], walk=9,
                                                                           step=1)
            else:
                self.id_grepseq_dic[f'{f_id}BCprimer_R'] = None
            if uniseq:
                steps = len(uniseq) // 20
                self.id_grepseq_dic[f'{f_id}uniseq'] = self.mk_grepseq(ref_seq=uniseq, walk=17, step=steps)
            else:
                self.id_grepseq_dic[f'{f_id}uniseq'] = None

    def is_complete(self, f_id, read, read_rev_com, match):
        gerpseqs_leftseq = self.id_grepseq_dic[f'{f_id}leftseq']
        gerpseqs_rightseq = self.id_grepseq_dic[f'{f_id}rightseq']
        grepseqs_BCprimerF = self.id_grepseq_dic[f'{f_id}BCprimer_F']
        grepseqs_BCprimerR = self.id_grepseq_dic[f'{f_id}BCprimer_R']
        grepseqs_uniseq = self.id_grepseq_dic[f'{f_id}uniseq']

        left_range = 150
        right_range = 150

        if self.l_seq:
            left_range = len(self.l_seq)
        
        if self.r_seq:
            right_range = len(self.r_seq)

        match_nf_leftseq = self.osearch(seqs_used=gerpseqs_leftseq, seq_searched=read[1][:left_range], match=match)
        match_nr_leftseq = self.osearch(seqs_used=gerpseqs_leftseq, seq_searched=read_rev_com[1][:left_range], match=match)
        if match_nf_leftseq:
            match_nf_rightseq = self.osearch(seqs_used=gerpseqs_rightseq, seq_searched=read[1][-right_range:], match=match)
            if match_nf_rightseq:
                match_nf_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read[1][:21],
                                                  match=match)
                if match_nf_BCprimerF:
                    match_nf_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR, seq_searched=read[1][-20:],
                                                      match=match)
                    if match_nf_BCprimerR:
                        match_nf_uniseq = self.osearch(seqs_used=grepseqs_uniseq, seq_searched=read[1][left_range:-right_range],
                                                       match=match)
                        if match_nf_uniseq:
                            return "NF"
        elif match_nr_leftseq:
            match_nr_rightseq = self.osearch(seqs_used=gerpseqs_rightseq, seq_searched=read_rev_com[1][-right_range:],
                                             match=match)
            if match_nr_rightseq:
                match_nr_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read_rev_com[1][:21],
                                                  match=match)
                if match_nr_BCprimerF:
                    match_nr_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR,
                                                      seq_searched=read_rev_com[1][-20:],
                                                      match=match)
                    if match_nr_BCprimerR:
                        match_nr_uniseq = self.osearch(seqs_used=grepseqs_uniseq,
                                                       seq_searched=read_rev_com[1][left_range:-right_range],
                                                       match=match)
                        if match_nr_uniseq:
                            return "NR"

    def demultiplex(self, fastq_file, match=1):
        logger.info(f'Demultiplexing {fastq_file}')
        nt = 0
        start = time.time()
        for read in read_fastq(fastq_file):
            # Obtain the reverse complementary sequence
            revcom_line = self.seq_rev_com(read)
            for f_id in self.fastq_ids:
                result = self.is_complete(f_id, read, revcom_line, match)

                if result == "NF":
                    line = f"{read[0]}\n{read[1]}\n{read[2]}\n{read[3]}\n"
                    self.id_reads_dic[f_id].append(line)

                elif result == "NR":
                    line_rc = f"{revcom_line[0]}\n{revcom_line[1]}\n{revcom_line[2]}\n{revcom_line[3]}\n"
                    self.id_reads_dic[f_id].append(line_rc)

            nt += 1
            if nt % 100000 == 0:
                logger.info("Processed %d reads in %.1f minutes" % (nt, (time.time() - start) / 60))

    def write_fastq(self):
        self.demulti_outpath = self.set_dir('Demultiplexed')
        self.demultis, stats = [], []
        self.demulti_random200s = {}
        stats = ['file\treads\n']
        for id, reads in self.id_reads_dic.items():
            out = os.path.join(self.demulti_outpath, id + '.fastq')
            out_random200 = os.path.join(self.demulti_outpath, id + '_random200.fastq')
            self.demultis.append(out)
            self.demulti_random200s[out] = None
            with open(out, 'w') as f_out:
                if reads:
                    lines = reads[:-2]
                    lines.append(reads[-1].strip())
                    f_out.writelines(lines)
                    if len(lines) > 200:
                        self.demulti_random200s[out] = out_random200
                        random_lines = random.sample(lines, 200)
                        with open(out_random200, 'w') as random_out:
                            random_out.writelines(random_lines)
                else:
                    continue
            logger.info(f"{len(reads)} reads are written to the {out}")
            stats.append(f"{out}\t{len(reads)}\n")
            seqkit_watch_cmd = f'seqkit watch -Q --fields ReadLen {out} -O {out.strip(".fastq")}.pdf'
            subprocess.call(seqkit_watch_cmd, shell=True)
        stats_out = os.path.join(self.demulti_outpath, 'Demultiplex_stats.txt')
        with open(stats_out, 'w') as f_stat:
            f_stat.writelines(stats)

    def visualization(self):
        self.visua_abspath = self.set_dir('Visualization')
        sorted_bams = []
        for demulti_id, demulti in zip(self.fastq_ids, self.demultis):
            ref_id = self.id_ref[demulti_id]
            reference = self.id_references[ref_id]
            demulti_random200 = self.demulti_random200s[demulti]

            sorted_bam = vis.visualizing(reference, demulti, self.visua_abspath, demulti_random200)
            sorted_bams.append(sorted_bam)
        logger.info("Visualization done")


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')
    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('-d', '--demultiplexinfo', help="Specify the DemultiplexInfo file", required=True)
    all_parser.add_argument('-n', '--nanopore', help="Specify the Nanopore sequencing data file", required=True)
    all_parser.add_argument('-r', '--referenceinfo', help="Specify the ReferenceInfo file", required=True)

    preprocess_parser = subparsers.add_parser('mkreference',
                                              help='Create FASTA files based on the information in the "ReferenceInfo"')
    preprocess_parser.add_argument('-r', '--referenceinfo', help="Specify the ReferenceInfo file", required=True)

    demultiplex_parser = subparsers.add_parser('demultiplex',
                                               help='Demultiplex FASTQ files based on the information in the "DemultiplexInfo"')
    demultiplex_parser.add_argument('-n', '--nanopore', help="Specify the Nanopore sequencing data file", required=True)
    demultiplex_parser.add_argument('-d', '--demultiplexinfo', help="Specify the DemultiplexInfo file", required=True)

    visualize_parser = subparsers.add_parser('visualize', help='Minimap2 align FASTQ files, '
                                                               'Samtools generate the sorted.bam and bai files '
                                                               'that IGV needed')
    visualize_parser.add_argument('-a', '--fasta', help='Specify the reference FASTA file', required=True)
    visualize_parser.add_argument('-q', '--fastq', help='Specify demultiplexed FASTQ files', required=True, nargs='+')

    return parser.parse_args()


def main():
    def if_list(var):
        if isinstance(var, list):
            return var
        else:
            return [var]

    args = parse_args()
    if args.command == 'all':
        g = GREPoreSeq()
        g.input_info(ref_yaml=args.referenceinfo, demulti_yaml=args.demultiplexinfo)
        g.mk_reference()
        g.disassemble()
        g.prepare_grepseq()
        g.demultiplex(fastq_file=args.nanopore, match=2)
        g.write_fastq()
        g.visualization()


    elif args.command == 'mkreference':
        """
        Run just mkreference step given the RefInfo
        """
        g = GREPoreSeq()
        g.input_info(ref_yaml=args.referenceinfo)
        g.mk_reference()


    elif args.command == 'demultiplex':
        g = GREPoreSeq()
        g.input_info(demulti_yaml=args.demultiplexinfo)
        g.disassemble()
        g.prepare_grepseq()
        g.demultiplex(fastq_file=args.nanopore, match=2)
        g.write_fastq()

    elif args.command == 'visualize':
        g = GREPoreSeq()
        outpath = g.set_dir('Visualization')
        fastqs = if_list(args.fastq)
        for demuti in fastqs:
            vis.visualizing(args.fasta, demuti, outpath)


if __name__ == '__main__':
    main()
