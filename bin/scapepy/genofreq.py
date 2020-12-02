#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 13:49:19 2020

@author: Fred Jaya
"""

import re
from collections import Counter

def match_gt(gt):
    """
    Regex for genotype
    """
    r = re.compile("(\d+\/\d+)|(\.\/\.)|(\d+\|\d+)|(\.\|\.)")
    matched_gts = r.match(gt).group(0)
    # Phased and unphased genotypes treated the same
    matched_gts = re.sub("\/|\|", "_", matched_gts)
    return(matched_gts)

def count_genotypes(gt_info):
    """
    Get count of unique genotypes per sample
    """
    l = []
    
    for i in gt_info:
        sample_gt = match_gt(i)
        l.append(sample_gt)
    
    # Lazy test 
    #if len(l) != 15:
    #    print("ERROR: Missing genotype samples, ", len(l), "/15")
    #    exit()
    #print(l)
    return dict(Counter(l))

def site_pass(gt_counts):
    """
    Identify sites that have a single unique genotype
    """
    pass_filter = False
    
    for key, value in gt_counts.items():
        if value == 1 and not key == '._.':
            pass_filter = True
            
    return pass_filter

def write_site(writer, pass_filter, line):
    """
    Write site to new file if passes filter
    """
    if pass_filter:
        writer.writelines(line)
    return

def genofreq(gt_info, writer, line):
    """
    Retain only sites where a unique genotype appears once
    """
    gt_counts = count_genotypes(gt_info)
    pass_filter = site_pass(gt_counts)
    write_site(writer, pass_filter, line)
