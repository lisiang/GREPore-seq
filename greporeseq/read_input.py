# -*- coding: utf-8 _*_
# @Time    : 2021/12/20 15:50
# @Author  : SiangLi
# @File    : read_input.py
# @Desc    : Reading input information
import logging

import yaml

logger = logging.getLogger('root')
logger.propagate = False

def mk_greps(string, walk, step):
    s, length = 0, len(string)
    greps = [string[s + step * i:walk + step * i] for i in range(round(length / step))]
    grepws = [seq for seq in greps if len(seq) == walk]
    return grepws


def r_input(yaml_file):
    # Check file format
    if not yaml_file.endswith('.yaml'):
        logger.error('Please provide the file in the correct format')
        return
    with open(yaml_file, "r", encoding="utf8") as f:
        input_info = yaml.load(f, Loader=yaml.FullLoader)
        logger.info(f'Successful loading "{yaml_file}"')
    return input_info



