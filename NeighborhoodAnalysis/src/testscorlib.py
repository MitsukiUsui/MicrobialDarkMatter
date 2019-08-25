#!/usr/bin/env python3

import unittest

import numpy as np
from scorelib import score_naive, score_independent, score_conditional


class TestSegmentManager(unittest.TestCase):
    indicator_matrix = np.array([
        [0, 1, 0, 0, 0],
        [-1, 1, 0, 0, -1],
        [0, 1, 0, 0, -1],
        [0, 1, 0, -1, -1],
        [0, 0, 0, 0, 1],
        [-1, -1, 0, 0, 1],
        [0, 0, 0, 0, 0],
        [-1, 0, 0, 0, -1],
        [-1, -1, 0, 0, 0],
        [-1, -1, 0, -1, -1]
    ]).astype(float)

    def test_naive(self):
        self.assertEqual(
            score_naive(self.indicator_matrix),
            6 / 10
        )

    def test_independent(self):
        self.assertEqual(
            score_independent(self.indicator_matrix),
            1 - (1 - 0 / 5) * (1 - 4 / 7) * (1 - 0 / 10) * (1 - 0 / 8) * (1 - 2 / 5)
        )

    def test_conditional(self):
        self.assertEqual(
            score_conditional(self.indicator_matrix),
            (1 * 6 + 0 + 1 / 2 + 4 / 7 + 11 / 14) / 10
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
