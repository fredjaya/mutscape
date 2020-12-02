#!/usr/bin/env python3

import unittest
from genofreq import match_gt, site_pass
from orphan import get_worker_gt_info

class test_genofreq(unittest.TestCase):
    def setUp(self):
        """ match_gt """
        self.gt_normal = "0/1:3,2:5:75:75,0,120"
        self.gt_two_digits = "10/10:15:59"
        
        """ Site pass """
        # This site in .vcf was removed when it shouldn't have
        self.false_neg = {'0/1':14, '1/1':1}
        self.true_pos1 = {'./.':1, '1/1':14}
        self.true_pos2 = {'1/1':1, '0/0':1, '0/1':12}
        
        # Add phased genotypes
        self.false_neg_ph = {'0|1':14, '1/1':1}
        self.true_pos1_ph = {'.|.':1, '1/1':14}
        self.true_pos2_ph = {'1|1':1, '0/0':1, '0/1':12}
        
    def test_match_gt(self):
        self.assertEqual(
                match_gt(self.gt_normal), "0_1")

        self.assertEqual(
                match_gt(self.gt_two_digits), "10_10")
        
    def test_site_pass(self):
        self.assertTrue(
                site_pass(self.false_neg))
        
        self.assertTrue(
                site_pass(self.true_pos1))
        
        self.assertTrue(
                site_pass(self.true_pos2))
        
        self.assertTrue(
                site_pass(self.false_neg_ph))
        
        self.assertTrue(
                site_pass(self.true_pos1_ph))
        
        self.assertTrue(
                site_pass(self.true_pos2_ph))
        

class test_orphan(unittest.TestCase):
    def setUp(self):
        self.remove = ['./.:1,2:3:12:76,0,12', './.:.:.:.:.', './.:.:.:.:.', 
             './.:20:20:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.',
             './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.',
             './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.\n']
        self.retain =  ['./.:7:7:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.',
             '0/1:5,2:7:69:69,0,204', './.:.:.:.:.', './.:.:.:.:.',
             './.:.:.:.:.', './.:.:.:.:.', '0/1:6,2:8:66:66,0,246',
             '0/1:6,3:9:99:108,0,243', '0/1:7,2:9:63:63,0,288',
             '0/1:7,4:11:99:147,0,282', './.:.:.:.:.', '0/1:4,2:6:72:72,0,162\n']
        self.remove_ph = ['./.:1,2:3:12:76,0,12', './.:.:.:.:.', './.:.:.:.:.', 
             './.:20:20:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.',
             './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.',
             './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.', '.|.:.:.:.:.\n']
        self.retain_ph =  ['./.:7:7:.:.', './.:.:.:.:.', './.:.:.:.:.', './.:.:.:.:.',
             '0/1:5,2:7:69:69,0,204', './.:.:.:.:.', './.:.:.:.:.',
             './.:.:.:.:.', './.:.:.:.:.', '0/1:6,2:8:66:66,0,246',
             '0/1:6,3:9:99:108,0,243', '0/1:7,2:9:63:63,0,288',
             '0/1:7,4:11:99:147,0,282', './.:.:.:.:.', '0|1:4,2:6:72:72,0,162\n']

    def test_get_worker_gt_info(self):
        self.assertEqual(
                get_worker_gt_info(self.remove), ['./.:.:.:.:.\n'])

        self.assertEqual(
                get_worker_gt_info(self.retain), ['0/1:4,2:6:72:72,0,162\n'])

        self.assertEqual(
                get_worker_gt_info(self.remove_ph), ['.|.:.:.:.:.\n'])

        self.assertEqual(
                get_worker_gt_info(self.retain_ph), ['0|1:4,2:6:72:72,0,162\n'])

if __name__ == '__main__':
    unittest.main()

