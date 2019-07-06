#!/usr/bin/env python3

from collections import defaultdict

import pandas as pd
from scipy import stats


class SegmentManager:
    def __init__(self):
        self.next_id = 0
        self.segment2members = dict() #key: segment_id, val: list of members
        self.size2segments = defaultdict(set) #key: size, val: set of segment_ids of the size
        self.member_count = 0

    def __len__(self):
        return len(self.segment2members)

    def add(self, members):
        segment_id = self.next_id
        size = len(members)
        self.next_id += 1
        self.segment2members[segment_id] = members
        self.size2segments[size].add(segment_id)
        self.member_count += size
        return segment_id

    def delete(self, segment_id):
        size = len(self.segment2members[segment_id])
        del self.segment2members[segment_id]
        self.size2segments[size].remove(segment_id)
        self.member_count -= size

    def split(self, segment_id, idx):
        members = self.get_members_by_id(segment_id)
        self.delete(segment_id)
        segment_id1 = self.add(members[:idx])
        segment_id2 = self.add(members[idx:])
        return segment_id1, segment_id2

    def get_members_by_id(self, segment_id):
        return self.segment2members[segment_id]

    def get_segments_by_size(self, size):
        return self.size2segments[size]

    def get_max_segment_size(self):
        return max(self.size2segments.keys())

    def get_segment_count(self):
        return len(self.segment2members)

    def get_member_count(self):
        return self.member_count

    def to_wcf(self):
        x = []
        y = []
        for size, segment_ids in self.size2segments.items():
            if len(segment_ids) > 0:
                x.append(size)
                y.append(len(segment_ids))
        wcf = Wcf(x, y)
        return wcf

    def to_df(self):
        if len(self.segment2members) == 0:
            return pd.DataFrame(columns=["segment_id", "member"])

        records = []
        for segment_id, members in self.segment2members.items():
            records += [{"segment_id": segment_id, "member": m} for m in members]
        df = pd.DataFrame(records, columns=["segment_id", "member"])
        return df

class Wcf:
    """
    weighted cumlative frequencies
    """
    def __init__(self, x, y):
        assert len(x) == len(y)
        weights = [x*y for x, y in zip(x, y)]
        res = stats.cumfreq(x, numbins=max(x)+1, defaultreallimits=(0, max(x)+1), weights=weights)
        self.wcf = res.cumcount / res.cumcount[-1]

    def __getitem__(self, key):
        if key < len(self.wcf):
            return self.wcf[key]
        else:
            return 1.0

    def __len__(self):
        return len(self.wcf)

    def to_array(self):
        return self.wcf
