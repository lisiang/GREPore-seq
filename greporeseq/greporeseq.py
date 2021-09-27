# -*- coding: utf-8 -*-
# @Time         : 2021/7/27
# @Author       : Siang_Li
# @File         : greporeseq_ori.py
# @Software     : PyCharm
# @E-mail       : lisiang.darcy@gmail.com
# @Description  : a pipeline for analysing ONT sequencing data

import os
import argparse

import log
import preprocess as pre
import demultiplex as dem
import visualization as vis
import stats as sta
from cigar import TargetRatio

logger = log.createCustomLogger('root')


def rm_dup_none(l):
    return list(filter(None, sorted(set(l), key=l.index)))


def filter_none(l):
    return list(filter(None, l))


class GREPoreSeq:
    def __init__(self):
        pass

    def set_dir(self, wd, dir_name):
        """Setting directory"""
        dir_abspath = os.path.join(wd, dir_name)
        if dir_name not in os.listdir(wd):
            os.mkdir(dir_abspath)
        return dir_abspath

    def excel_info(self, excel_name):
        """Extraction of excel information"""
        logger.info("Loading {0}...".format(excel_name))
        self.site_names = pre.read_excel(excel_name, 0)
        self.site_seqs = pre.read_excel(excel_name, 1, upper=True)
        self.bc_lenths = list(filter(None, pre.read_excel(excel_name, 4, dup_rm=False)))
        self.bc_primer_names = pre.read_excel(excel_name, 5, dup_rm=False)
        self.bc_primer_seqs = pre.read_excel(excel_name, 6, dup_rm=False, upper=True)

    def consolidate_np(self, NP_data_name):
        """Forward and reverse consolidation of NP data"""
        chopped_fr_fastq_gz = pre.consolidate_np_data(NP_data_name)
        self.chopped_fr_name = chopped_fr_fastq_gz

    def make_bc_primer_seq(self):
        """according to excel info make BC-primer-seq file"""
        self.pre_abspath = self.set_dir(os.getcwd(), 'preprocess')
        try:
            self.bcps_files = pre.make_bcps(self.bc_primer_names, self.bc_primer_seqs, self.bc_lenths, self.pre_abspath)
        except Exception as e:
            logger.error(e)
        else:
            logger.info('BC-primer-seq file has been generated')

    def make_grepseq(self):
        """according to Amplicon fa file generate Grepseq file"""
        self.pre_abspath = self.set_dir(os.getcwd(), 'preprocess')

        try:
            self.amp_fa_list = pre.make_fa(self.site_names, self.site_seqs, self.pre_abspath)
        except Exception as e:
            logger.error(e)
        else:
            logger.info('AMP fa file has been generated')

        self.grep_seq_files = pre.make_grepseq(self.amp_fa_list)
        logger.info('Grepseq file has been generated')
        self.grep_seqL = sorted([name for name in self.grep_seq_files if 'GrepL' in name],
                                key=self.grep_seq_files.index)
        self.grep_seqR = sorted([name for name in self.grep_seq_files if 'GrepR' in name],
                                key=self.grep_seq_files.index)

    def isolate(self):
        """demultiplex the NP-FR data"""
        self.dem_abspath = self.set_dir(os.getcwd(), 'demultiplexed')
        try:
            self.pcr_products = dem.isolate_PCR_products(self.chopped_fr_name,
                                                         self.grep_seqL,
                                                         self.grep_seqR,
                                                         self.dem_abspath)
        except Exception as e:
            logger.error(e)

    def demultiplex(self):
        """demultiplex the isolated datat"""
        self.dem_abspath = self.set_dir(os.getcwd(), 'demultiplexed')
        try:
            demulted_info = dem.demultiplex_PCR_products(self.pcr_products,
                                                         self.bcps_files,
                                                         self.dem_abspath)
            self.demulted_products = demulted_info[0]
            self.site_demulted_dict = demulted_info[1]
        except Exception as e:
            logger.error(e)
        if self.demulted_products:
            logger.info('Demultiplex has been done')
        else:
            logger.error("Please ensure {}'s name contains  {}'s info ".format(self.bcps_files,
                                                                               self.pcr_products))

    def visualization(self):
        self.visua_abspath = self.set_dir(os.getcwd(), 'visualization')
        try:
            self.mmi_file_list = vis.make_mmi_file(self.amp_fa_list,
                                                   self.visua_abspath)
            self.sorted_bam_list = vis.visualizing(self.mmi_file_list,
                                                   self.demulted_products,
                                                   self.visua_abspath)

        except Exception as e:
            logger.error(e)
        if self.sorted_bam_list:
            pass
        else:
            logger.error("Please ensure fq's name contain fa file's info ")

    def stats_all_fq(self):
        """stat all fastq file's read info"""
        self.dem_abspath = self.set_dir(os.getcwd(), 'demultiplexed')
        self.stat_abspath = self.set_dir(os.getcwd(), 'statistical-results')
        try:
            sta.stats(self.dem_abspath, self.stat_abspath)
        except Exception as e:
            logger.error(e)
        else:
            logger.info('Stats done')

    def conculate_coverage(self):
        """calculate the files coverage"""
        self.visua_abspath = self.set_dir(os.getcwd(), 'visualization')
        self.stat_abspath = self.set_dir(os.getcwd(), 'statistical-results')
        try:
            sta.coverage(self.visua_abspath, self.stat_abspath)
        except Exception as e:
            logger.error(e)
        else:
            logger.info('Calculation done')

    def stat_cigar(self, *sorted_bam, target, length, outname):
        """stat cigar string info"""
        self.stat_abspath = self.set_dir(os.getcwd(), 'statistical-results')
        output = os.path.join(self.stat_abspath, outname)
        t = TargetRatio()
        t.write_ratio(*sorted_bam, target=target, length=length, output=output)

    def large_deletion(self, *bam):
        """stat CIGAR-D 100 analyse large deletion ratio"""
        sorted_bam = []
        if bam:
            sorted_bam = bam
        else:
            sorted_bam = self.sorted_bam_list
        logger.info('Calculating large deletion')
        try:
            self.stat_cigar(*sorted_bam, target='D', length=100, outname='Large-deletion-ratio.txt')
        except Exception as e:
            logger.error(e)
        else:
            logger.info('Calculation done')


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(description='Individual Step Commands',
                                       help='Use this to run individual steps of the pipeline',
                                       dest='command')
    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline')
    all_parser.add_argument('-e', '--excel', help='Specify the excel name', required=True)
    all_parser.add_argument('-N', '--Nanopore', help="Specify the Nanopore sequencing data's name",
                            required=True)
    all_parser.add_argument('-l', '--large', help="Calculate the Large deletion ratio of demultiplexed files",
                            action='store_true')
    all_parser.add_argument('-c', '--coverage', help="Use samtools to calculate coverage",
                            action='store_true')

    preprocess_parser = subparsers.add_parser('preprocess', help='Perform preprocessing steps')
    preprocess_parser.add_argument('-e', '--excel', help='Specify the excel name', required=True)
    preprocess_parser.add_argument('-N', '--Nanopore', help="Specify the Nanopore sequencing data's name",
                                   required=True)

    isolate_parser = subparsers.add_parser('isolate', help='Using Grepseq isolate PCR products')
    isolate_parser.add_argument('-N', '--NanoporeData', help='Nanopore datat after preprocessing', required=True)
    isolate_parser.add_argument('-lg', '--left_grep', help='Specify the Grepseq left file', required=True)
    isolate_parser.add_argument('-rg', '--right_grep', help='Specify the Grepseq right file', required=True)

    demultiplex_parser = subparsers.add_parser('demultiplex', help='Demultiplex isolated FASTQ files')
    demultiplex_parser.add_argument('-p', '--products', help='Specify the isolated FASTQ files', required=True,
                                    nargs='+')
    demultiplex_parser.add_argument('-b', '--BC_primer_seq', help='Specify the BC-primer-seq file', required=True,
                                    nargs='+')

    visualize_parser = subparsers.add_parser('visualize', help='Minimap2 align FASTQ files, '
                                                               'Samtools generate the sorted.bam and bai files '
                                                               'that IGV needed')
    visualize_parser.add_argument('-a', '--fasta', help='Specify the fasta file', required=True, nargs='+')
    visualize_parser.add_argument('-i', '--input', help='Specify the demultiplexed file', required=True, nargs='+')

    stats_parser = subparsers.add_parser('stats', help='Stats preprocessed and demultiplexed FASTQ files')

    coverage_parser = subparsers.add_parser('coverage', help='Use samtools to calculate coverage')

    large_del_parser = subparsers.add_parser('large', help='Calculate the Large deletion ratio of demultiplexed files')
    large_del_parser.add_argument('-b', '--bam',
                                  help='The bam format file that needs to calculate the large deletion ratio',
                                  nargs='+', required=True)

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
        g.excel_info(args.excel)
        g.consolidate_np(args.Nanopore)
        g.make_bc_primer_seq()
        g.make_grepseq()
        g.isolate()
        g.demultiplex()
        g.visualization()
        g.stats_all_fq()
        if args.large:
            g.large_deletion()
        if args.coverage:
            g.conculate_coverage()


    elif args.command == 'preprocess':
        """
        Run just preprocessing step given the excel_file and Nanopore_Data
        """
        g = GREPoreSeq()
        g.excel_info(args.excel)
        g.consolidate_np(args.Nanopore)
        g.make_bc_primer_seq()
        g.make_grepseq()
        if args.large:
            g.large_deletion()

    elif args.command == 'isolate':
        g = GREPoreSeq()
        g.chopped_fr_name = args.NanoporeData
        g.grep_seqL = [args.left_grep]
        g.grep_seqR = [args.right_grep]
        g.isolate()

    elif args.command == 'demultiplex':
        g = GREPoreSeq()
        g.pcr_products = if_list(args.bam)
        g.bcps_files = if_list(args.BC_primer_seq)
        g.demultiplex()

    elif args.command == 'visualize':
        g = GREPoreSeq()
        g.amp_fa_list = if_list(args.fasta)
        g.demulted_products = if_list(args.input)
        g.visualization()

    elif args.command == 'stats':
        g = GREPoreSeq()
        g.stats_all_fq()

    elif args.command == 'coverage':
        g = GREPoreSeq()
        g.conculate_coverage()

    elif args.command == 'large':
        g = GREPoreSeq()
        bams = if_list(args.bam)
        g.large_deletion(*bams)


if __name__ == '__main__':
    main()
