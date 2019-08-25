#!/usr/bin/env python3

import unittest

from splitlib import SegmentManager, Wcf


class TestSegmentManager(unittest.TestCase):
    def test_add(self):
        segment_manager = SegmentManager()
        self.assertEqual(segment_manager.get_segment_count(), 0)
        self.assertEqual(segment_manager.get_member_count(), 0)

        member1 = [1, 2, 3]
        segment_id1 = segment_manager.add(member1)
        self.assertEqual(segment_manager.get_segment_count(), 1)
        self.assertEqual(segment_manager.get_member_count(), len(member1))
        self.assertEqual(segment_manager.get_members_by_id(segment_id1), member1)

        member2 = [5, 6]
        segment_id2 = segment_manager.add(member2)
        self.assertEqual(segment_manager.get_segment_count(), 2)
        self.assertEqual(segment_manager.get_member_count(), len(member1) + len(member2))
        self.assertEqual(segment_manager.get_members_by_id(segment_id1), member1)
        self.assertEqual(segment_manager.get_members_by_id(segment_id2), member2)

    def test_delete(self):
        segment_manager = SegmentManager()
        member1 = [1, 2, 3]
        segment_id1 = segment_manager.add(member1)
        member2 = [5, 6]
        segment_id2 = segment_manager.add(member2)
        segment_manager.delete(segment_id1)
        self.assertEqual(segment_manager.get_segment_count(), 1)
        self.assertEqual(segment_manager.get_member_count(), len(member2))
        self.assertEqual(segment_manager.get_members_by_id(segment_id2), member2)

    def test_split(self):
        segment_manager = SegmentManager()
        member = [1, 2, 3, 4, 5]
        segment_id = segment_manager.add(member)
        idx = 2
        child_id1, child_id2 = segment_manager.split(segment_id, idx)
        self.assertEqual(segment_manager.get_segment_count(), 2)
        self.assertEqual(segment_manager.get_member_count(), len(member))
        self.assertEqual(segment_manager.get_members_by_id(child_id1), member[:idx])
        self.assertEqual(segment_manager.get_members_by_id(child_id2), member[idx:])


class TestCwf(unittest.TestCase):
    def test_cwf(self):
        x = [1, 2, 4]
        y = [4, 1, 1]
        wcf = Wcf(x, y)
        self.assertEqual(wcf[0], 0.0)
        self.assertEqual(wcf[1], 0.4)
        self.assertEqual(wcf[2], 0.6)
        self.assertEqual(wcf[3], 0.6)
        self.assertEqual(wcf[4], 1.0)
        self.assertEqual(wcf[5], 1.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
