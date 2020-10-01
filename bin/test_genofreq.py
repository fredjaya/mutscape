#!/usr/bin/env python3

import unittest
from genofreq import match_gt
from genofreq import get_worker_gt_info

class test_genofreq(unittest.TestCase):
    def setUp(self):
        self.gt_normal = "0/1:3,2:5:75:75,0,120"
        self.gt_two_digits = "10/10:15:59"

    def test_match_gt(self):
        self.assertEqual(
                match_gt(self.gt_normal),
                "0/1")

        self.assertEqual(
                match_gt(self.gt_two_digits),
                "10/10")

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

    def test_get_worker_gt_info(self):
        self.assertTrue(
                get_worker_gt_info(self.remove),
                ['./.:.:.:.:.\n'])

        self.assertTrue(
                get_worker_gt_info(self.retain),
                ['0/1:4,2:6:72:72,0,162\n'])

if __name__ == '__main__':
    unittest.main()

