# -*- coding: utf-8 -*-
"""一个在命令窗口显示信息的模块"""
import logging


def createCustomLogger(name):
    formatter = logging.Formatter(fmt='[%(asctime)s][%(levelname)s][%(module)s] %(message)s',
                                  datefmt='%m/%d %I:%M:%S%p')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    return logger
