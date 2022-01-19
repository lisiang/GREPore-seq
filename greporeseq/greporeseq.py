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
        self.fastq_ids, self.left_150s, self.right_150s = [], [], []
        self.id_ref, self.BCprimer_Fs, self.BClen_Fs, self.id_BCprimer_Rs, self.id_uniseq, self.id_BClen_R = {}, {}, {}, {}, {}, {}
        for demuti_id, demuti_info in self.demul_info.items():
            demuti_id = self.reform_id(demuti_id)
            self.fastq_ids.append(demuti_id)
            self.left_150s.append(self.reform_seq(demuti_info['left_150bp']))
            self.right_150s.append(self.reform_seq(demuti_info['right_150bp']))
            self.BCprimer_Fs[f'{demuti_id}'] = self.reform_seq(demuti_info['BCprimer_F'])
            self.BClen_Fs[f'{demuti_id}'] = self.reform_num(demuti_info['BClen_F'])
            self.id_ref[f'{demuti_id}'] = self.reform_id(demuti_info['reference_id'])
            self.id_BCprimer_Rs[f'{demuti_id}'] = self.rev_com(self.reform_seq(demuti_info['BCprimer_R']))
            self.id_BClen_R[f'{demuti_id}'] = self.reform_num(demuti_info['BClen_R'])
            self.id_uniseq[f'{demuti_id}'] = self.reform_seq(demuti_info['unique_sequence'])

    def prepare_grepseq(self):
        logger.info('Preparing Grepseqs')
        self.id_reads_dic, self.id_grepseq_dic = {}, {}
        for f_id, l_150, r_150 in zip(self.fastq_ids, self.left_150s, self.right_150s):
            BCprimer_F = self.BCprimer_Fs[f_id]
            BCprimerR = self.id_BCprimer_Rs[f_id]
            BClen_F = self.BClen_Fs[f_id]
            BClen_R = self.id_BClen_R[f_id]
            uniseq = self.id_uniseq[f_id]

            if not BClen_F:
                BCrange_F = -1
            else:
                BCrange_F = BClen_F + 8

            if not BClen_R:
                BCrange_R = -1
            else:
                BCrange_R = BClen_R + 8

            self.id_reads_dic[f_id] = []
            self.id_grepseq_dic[f'{f_id}left150'] = self.mk_grepseq(ref_seq=l_150[:151], walk=15, step=20)
            self.id_grepseq_dic[f'{f_id}right150'] = self.mk_grepseq(ref_seq=r_150[-150:], walk=15, step=20)

            if BCprimer_F:
                self.id_grepseq_dic[f'{f_id}BCprimer_F'] = self.mk_grepseq(ref_seq=BCprimer_F[:BCrange_F], walk=9,
                                                                           step=1)
            else:
                self.id_grepseq_dic[f'{f_id}BCprimer_F'] = None
            if BCprimerR:
                self.id_grepseq_dic[f'{f_id}BCprimer_R'] = self.mk_grepseq(ref_seq=BCprimerR[:BCrange_R], walk=9,
                                                                           step=1)
            else:
                self.id_grepseq_dic[f'{f_id}BCprimer_R'] = None
            if uniseq:
                self.id_grepseq_dic[f'{f_id}uniseq'] = self.mk_grepseq(ref_seq=uniseq, walk=15, step=200)
            else:
                self.id_grepseq_dic[f'{f_id}uniseq'] = None

    def is_complete(self, f_id, read, read_rev_com, match):
        gerpseqs_left150 = self.id_grepseq_dic[f'{f_id}left150']
        gerpseqs_right150 = self.id_grepseq_dic[f'{f_id}right150']
        grepseqs_BCprimerF = self.id_grepseq_dic[f'{f_id}BCprimer_F']
        grepseqs_BCprimerR = self.id_grepseq_dic[f'{f_id}BCprimer_R']
        grepseqs_uniseq = self.id_grepseq_dic[f'{f_id}uniseq']

        match_nf_left150 = self.osearch(seqs_used=gerpseqs_left150, seq_searched=read[1][:151], match=match)
        match_nr_left150 = self.osearch(seqs_used=gerpseqs_left150, seq_searched=read_rev_com[1][:151], match=match)
        if match_nf_left150:
            match_nf_right150 = self.osearch(seqs_used=gerpseqs_right150, seq_searched=read[1][-150:], match=match)
            if match_nf_right150:
                match_nf_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read[1][:21],
                                                  match=match)
                if match_nf_BCprimerF:
                    match_nf_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR, seq_searched=read[1][-20:],
                                                      match=match)
                    if match_nf_BCprimerR:
                        match_nf_uniseq = self.osearch(seqs_used=grepseqs_uniseq, seq_searched=read[1][151:-150],
                                                       match=match)
                        if match_nf_uniseq:
                            return "NF"
        elif match_nr_left150:
            match_nr_right150 = self.osearch(seqs_used=gerpseqs_right150, seq_searched=read_rev_com[1][-150:],
                                             match=match)
            if match_nr_right150:
                match_nr_BCprimerF = self.osearch(seqs_used=grepseqs_BCprimerF, seq_searched=read_rev_com[1][:21],
                                                  match=match)
                if match_nr_BCprimerF:
                    match_nr_BCprimerR = self.osearch(seqs_used=grepseqs_BCprimerR,
                                                      seq_searched=read_rev_com[1][-20:],
                                                      match=match)
                    if match_nr_BCprimerR:
                        match_nr_uniseq = self.osearch(seqs_used=grepseqs_uniseq,
                                                       seq_searched=read_rev_com[1][151:-150],
                                                       match=match)
                        if match_nr_uniseq:
                            return "NR"

    def demultiplex(self, fastq_file, match=2):
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
        stats = ['file\treads\n']
        for id, reads in self.id_reads_dic.items():
            out = os.path.join(self.demulti_outpath, id + '.fastq')
            self.demultis.append(out)
            with open(out, 'w') as f_out:
                if reads:
                    lines = reads[:-2]
                    lines.append(reads[-1].strip())
                    f_out.writelines(lines)
                else:
                    continue
            logger.info(f"{len(reads)} reads are written to the {out}")
            stats.append(f"{out}\t{len(reads)}\n")
        stats_out = os.path.join(self.demulti_outpath, 'Demultiplex_stats.txt')
        with open(stats_out, 'w') as f_stat:
            f_stat.writelines(stats)

    def visualization(self):
        self.visua_abspath = self.set_dir('Visualization')
        sorted_bams = []
        for demulti_id, demulti in zip(self.fastq_ids, self.demultis):
            ref_id = self.id_ref[demulti_id]
            reference = self.id_references[ref_id]

            sorted_bam = vis.visualizing(reference, demulti, self.visua_abspath)
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

    preprocess_parser = subparsers.add_parser('reference',
                                              help='Create FASTA files based on the information in the "RefenrenceInfo"')
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
