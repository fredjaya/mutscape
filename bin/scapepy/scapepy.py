#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:30 2020

@author: Fred Jaya
"""

import argparse
from genofreq import genofreq
from orphan import remove_orphan_sites
from single import remove_single_sites
from uniqueworker import remove_unique_workers

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

                    elif apply_filter == "single":
                        remove_single_sites(gt_info, writer, line)

                    elif apply_filter == "uniqueworker":
                        remove_unique_workers(gt_info, writer, line)
    return

if __name__ == "__main__":

    """
    Set arguments
    """                                                       
    parser = argparse.ArgumentParser()                        
    parser.add_argument("filter", help = "genofreq | orphan | single | uniqueworker")
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
    elif args.filter == "single":
        print("Removing sites where SNP is called in one sample only")
    elif args.filter == "uniqueworker":
        print("Removing sites where Worker genotype is unique")
                                                                          
    parse_vcf(args.filter, args.vcf_input, args.vcf_output)   
    
