#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 13:49:19 2020

@author: Fred Jaya
"""

import argparse
import os
import re
from collections import Counter

def count_genotypes(gt_info):
    """
    Get count of unique genotypes per sample
    """
    l = []
    
    for i in gt_info:
        sample_gt = match_gt(i)
        l.append(sample_gt)
    
    # Lazy test 
    if len(l) != 15:
        print("ERROR: Missing genotype samples, ", len(l), "/15")
        exit()

    return dict(Counter(l))

def match_gt(gt):
    """
    Regex for genotype
    """
    r = re.compile("(\d+\/\d+)|(\.\/\.)")
    return r.match(gt).group(0)

def site_pass(gt_counts):
    """
    Identify sites that have a single unique genotype
    """
    for key, value in gt_counts.items():
        if value == 1 and not key == "./.":
            return True
        return False 

def write_site(writer, pass_filter, line):
    """
    Write site to new file if passes filter
    """
    if pass_filter:
        writer.writelines(line)
    return

def parse_vcf(input_vcf, output_vcf):
    """
    Iterate through each variant and write to new file if passes filters
    """
    with open(input_vcf, 'r') as reader:
        with open(output_vcf, 'w') as writer:
            for line in reader:
                if line[0] == "#":
                    writer.writelines(line)
                else:
                    gt_info = line.split("\t")[9:]
                    gt_counts = count_genotypes(gt_info)
                    pass_filter = site_pass(gt_counts)
                    write_site(writer, pass_filter, line)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_input")
    parser.add_argument("vcf_output")
    args = parser.parse_args()
    parse_vcf(args.vcf_input, args.vcf_output)
