# -*- coding: utf-8 -*-
# @Time         : 2021/8/5
# @Author       : Siang_Li
# @File         : cigar.py
# @Software     : PyCharm
# @E-mail       : lisiang.darcy@gmail.com
# @Description  : stat CIGAR string info in sam format file
import os
import re

from tqdm import tqdm
import pysam


class TargetRatio:
    def __init__(self):
        self.name_dict = {}
        self.cigar_dict = {}
        self.name_cigar_dict = {}
        self.reduplicate_dict = {}
        self.specified_dict = {}
        self.filtered_num = 1

    def filter(self, bam_file):
        """filter the read which mapping quality lower than 30"""
        filtered_num = 1
        wb = pysam.AlignmentFile(bam_file, 'rb')
        for read in wb:
            if read.mapping_quality != 255 and read.mapping_quality >= 30:
                filtered_num += 1
                yield read
        self.filtered_num = filtered_num
        wb.close()

    def get_bam_info(self, bam):
        """read bam file information"""
        read_names = []
        read_cigars = []
        for read in self.filter(bam):
            read_names.append(list(str(read).split('\t'))[0])
            read_cigars.append(list(str(read).split('\t'))[5])
        series = range(1, self.filtered_num)
        name_dict = dict(zip(series, read_names))
        cigar_dict = dict(zip(series, read_cigars))
        return name_dict, cigar_dict

    def get_keys(self, d, value):
        """according to value take key"""
        return [k for k, v in d.items() if v == value]

    def collapse(self, dict):
        """collapse reduplicate name and specified name"""
        name_list = list(dict.values())
        set_list = set(name_list)
        reduplicate_dict = {}
        specified_dict = {}
        for name in set_list:
            if name_list.count(name) > 1:
                reduplicate_dict[name] = self.get_keys(dict, name)
            elif name_list.count(name) == 1:
                specified_dict[name] = self.get_keys(dict, name)
        return reduplicate_dict, specified_dict

    def combine_cigar(self, d, cigar_dict):
        """combine cigar string"""
        name_cigar = {}
        for name, rows in d.items():
            cigar_combine = ''
            for row in rows:
                cigar_combine += cigar_dict[row]
            name_cigar[name] = cigar_combine
        return name_cigar

    def write_ratio(self, *bam_files, target, length, output, Opposing=False):
        """write result to txt file"""
        with open(f'{output}.txt', 'w') as deletion_object:
            deletion_object.write(
                'name\tfiltered-reads(MAPQ>=30)\trmdup-reads\t'
                f'{target}>={length}bp\t{target}{length}ratio'
                '\n')
            inters = tqdm(bam_files)
            for bam_file in inters:
                bam_info = self.get_bam_info(bam_file)
                self.name_dict = bam_info[0]
                self.cigar_dict = bam_info[1]

                collapse_result = self.collapse(self.name_dict)
                self.reduplicate_dict = collapse_result[0]
                self.specified_dict = collapse_result[1]

                self.name_cigar_dict = (self.combine_cigar(self.reduplicate_dict, self.cigar_dict))
                self.name_cigar_dict.update(self.combine_cigar(self.specified_dict, self.cigar_dict))
                if not self.name_cigar_dict:
                    continue

                stat = 0
                for name, cigar_string in self.name_cigar_dict.items():
                    mode = "\\d{2,}"
                    re_string = mode + f'[{target}]'
                    Conformity = re.findall(re_string, cigar_string)
                    if Conformity:
                        numbers = [int(i[:i.index(f'{target}')]) for i in Conformity]
                        logs = [num for num in numbers if num >= length]
                        if logs:
                            stat += 1

                line = f"{bam_file}\t" \
                       f"{self.filtered_num}\t" \
                       f"{len(self.name_cigar_dict)}\t" \
                       f"{stat}\t"
                if Opposing:
                    line += f"{(1-(stat/len(self.name_cigar_dict)))*100:.2f}%\n"
                else:
                    line += f"{(stat / len(self.name_cigar_dict) * 100):.2f}%\n"
                deletion_object.write(line)
