#!/usr/bin/env python3

import unittest

class test_minGeno(unittest.TestCase):
    def setUp(self):
        self.gt_two_digits = "10/10:15:59"

    def test_match_gt(self):
        self.assertEqual(
                match_gt(self.gt_two_digits),
                "10/10")

if __name__ == '__main__':
    unittest.main()
