import numpy as np


def score_naive(indicator_matrix):
    total = indicator_matrix.shape[0]
    found = ((indicator_matrix == 1).sum(axis=1) > 0).sum()
    if total > 0:
        return found / total
    else:
        return 0


def score_independent(indicator_matrix):
    total_arr = (indicator_matrix >= 0).sum(axis=0)
    found_arr = (indicator_matrix == 1).sum(axis=0)
    msk = total_arr > 0
    if msk.sum() == 0:
        return 0
    else:
        freq_arr = found_arr[msk] / total_arr[msk]
        score = 1 - np.prod(1 - freq_arr)
        return score


def score_conditional(indicator_matrix):
    def calculate_phat(target_arr, score_arr):
        num = np.dot(target_arr == 1, 1 - score_arr)
        den = np.dot(target_arr >= 0, 1 - score_arr)
        if den > 0:
            return num / den
        elif den == 0:
            return 0

    score_arr = np.zeros(indicator_matrix.shape[0])
    total_arr = (indicator_matrix >= 0).sum(axis=0)
    for j in np.argsort(total_arr)[::-1]:
        target_arr = indicator_matrix[:, j].copy()
        phat = calculate_phat(target_arr, score_arr)
        target_arr[target_arr == -1] = phat  # fill missing value with estimated phat
        score_arr += target_arr * (1 - score_arr)
    score = np.mean(score_arr)
    return score
