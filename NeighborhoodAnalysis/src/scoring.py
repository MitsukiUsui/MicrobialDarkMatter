#!/usr/bin/env python3

import numpy as np


def score_naive(obs_mat, positions):
    total = obs_mat.shape[0]
    found = len(set([pos.i for pos in positions]))
    if total > 0:
        return found / total
    else:
        return 0


def score_closest(obs_mat, positions, calibration=True):

    def filter_closest(positions):
        closest = {} #key: i, val: pos
        for pos in positions:
            if pos.i not in closest or pos < closest[pos.i]:
                closest[pos.i] = pos
        return closest.values()

    score_mat = obs_mat.astype(int) - 1
    for pos in filter_closest(positions):
        score_mat[pos.i] = 0
        score_mat[pos.i][pos.j] = 1

    n0 = (score_mat==0).sum(axis=0) if calibration else (score_mat<=-0).sum(axis=0)
    n1 = (score_mat==1).sum(axis=0)
    return (n1 / (n0+n1)).sum()


def score_independent(obs_mat, positions, calibration=True):
    """
    Not equal to naive even if calibration = False
    """

    if calibration:
        trial_arr = obs_mat.sum(axis=0)
    else:
        trial_arr = obs_mat.shape[0] * np.ones(obs_mat.shape[1])

    count_arr = np.zeros(len(trial_arr))
    for pos in positions:
        count_arr[pos.j] += 1

    msk = trial_arr > 0
    freq_arr = count_arr[msk] / trial_arr[msk]
    score = 1 - np.prod(1 - freq_arr)
    return score


def score_conditional(obs_mat, positions, calibration=True):

    def caliculate_phat(target_vec, score_vec):
        num = np.dot(target_vec==1, 1-score_vec)
        den = np.dot(target_vec>=0, 1-score_vec)
        assert num >= 0 and den >= 0
        if den > 0:
            return num / den
        elif den == 0:
            return 0

    def score_conditional_js(obs_mat, positions, js, calibration=True):
        """
        caribrate with conditional probabilities in a given order of js
        """

        target_mat = obs_mat.astype(float) - 1
        for pos in positions:
            target_mat[pos.i][pos.j] = 1

        H, W = target_mat.shape
        score_vec = np.zeros(H)
        for j in js:
            target_vec = target_mat[:, j].copy()
            phat = caliculate_phat(target_vec, score_vec) if calibration else 0
            target_vec[target_vec==-1] = phat #fill missing value with estimated phat
            score_vec += target_vec * (1 - score_vec)

        score = np.mean(score_vec)
        return score

    H, W = obs_mat.shape
    assert W % 2 == 1
    DIST = int((W-1) / 2)
    js = list(np.argsort(obs_mat.sum(axis=0))[::-1])
    js.remove(DIST) #remove center position
    score = score_conditional_js(obs_mat, positions, js, calibration)
    return score

