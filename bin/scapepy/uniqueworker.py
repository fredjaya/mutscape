#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 16:42 2020

@author: Fred Jaya
"""

from genofreq import count_genotypes, write_site, match_gt
from orphan import get_worker_gt_info
import pandas as pd

def unique_worker(larv_counts, worker_gt):
    """
    Identify sites where Worker genotype is unique
    """
    pass_filter = False

    for key, value in larv_counts.items():
        if key == worker_gt:
            pass_filter = True

    return pass_filter

def remove_unique_workers(gt_info, writer, line):
    """
    Remove sites where Worker genotype is unique
    """
    larv_gt = gt_info[0:14]
    larv_counts = count_genotypes(larv_gt)
    
    worker_gt = get_worker_gt_info(gt_info)
    worker_gt = match_gt(worker_gt[0])

    pass_filter = unique_worker(larv_counts, worker_gt)

    write_site(writer, pass_filter, line)
    return
