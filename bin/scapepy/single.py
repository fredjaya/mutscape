#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:40 2020

@author: Fred Jaya
"""

from genofreq import count_genotypes, write_site

def is_single_site(gt_counts):
    """
    Identify sites that only have one site
    """
    pass_filter = True

    for key, value in gt_counts.items():
        if key == '._.' and value == 14:
            pass_filter = False

    return pass_filter

def remove_single_sites(gt_info, writer, line):
    """
    Remove sites where SNP is only called in 1 sample
    """
    gt_counts = count_genotypes(gt_info)
    pass_filter = is_single_site(gt_counts)
    write_site(writer, pass_filter, line)
    return
