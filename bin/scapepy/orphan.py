#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:40 2020

@author: Fred Jaya
"""

from genofreq import match_gt, write_site

def get_worker_gt_info(gt_info):
    return gt_info[14:15]

def retain_worker(gt):
    """
    Write site to new file if Worker genotype != ./.
    """
    if not gt == "._.":
        return True
    return False

def remove_orphan_sites(gt_info, writer, line):
    """
    Retain only sites where the Worker has a called a SNP
    """
    worker_gt_info = get_worker_gt_info(gt_info)
    worker_gt = match_gt(worker_gt_info[0])
    pass_filter = retain_worker(worker_gt)
    write_site(writer, pass_filter, line)
    return
