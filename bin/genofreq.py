#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 13:49:19 2020

@author: Fred Jaya
"""

import argparse
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
    print(matched_gts)
    return

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

def site_pass(gt_counts):
    """
    Identify sites that have a single unique genotype
    """
    pass_filter = False
    
    for key, value in gt_counts.items():
        if value == 1 and not key == "./.":
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

def get_worker_gt_info(gt_info):
    return gt_info[14:15]

def retain_worker(gt):
    """
    Write site to new file if Worker genotype != ./.
    """
    if not gt == "./.":
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

def parse_vcf(apply_filter, input_vcf, output_vcf):
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
                    
                    if apply_filter == "genofreq":
                        genofreq(gt_info, writer, line)
                    
                    elif apply_filter == "orphan":
                        remove_orphan_sites(gt_info, writer, line)
    return

if __name__ == '__main__':
    """
    Set arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("filter", help = "genofreq or orphan")
    parser.add_argument("vcf_input")
    parser.add_argument("vcf_output")
    args = parser.parse_args()
    
    """
    Run
    """
    if args.filter == "genofreq":
        print("Retaining sites with a unique genotype == 1")
    elif args.filter == "orphan":
        print("Retaining sites where Worker has a SNP called")

    parse_vcf(args.filter, args.vcf_input, args.vcf_output)
