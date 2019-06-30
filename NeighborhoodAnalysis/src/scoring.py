import numpy as np

def score_naive(indicator_matrix):
    total = indicator_matrix.shape[0]
    found = ((indicator_matrix==1).sum(axis=1) > 0).sum()
    if total > 0:
        return found / total
    else:
        return 0

def score_independent(indicator_matrix):
    total_arr = (indicator_matrix>=0).sum(axis=0)
    found_arr = (indicator_matrix==1).sum(axis=0)
    msk = total_arr > 0
    if msk.sum() == 0:
        return 0
    else:
        freq_arr = found_arr[msk] / total_arr[msk]
        score = 1 - np.prod(1 - freq_arr)
        return score

def score_conditional(gene_name, matrix):
    def indicator(position, gene_name):
        if position is None:
            return -1.0
        elif position.gene_name == gene_name:
            return 1.0
        else:
            return 0.0

    def caliculate_phat(target_vec, score_vec):
        num = np.dot(target_vec==1, 1-score_vec)
        den = np.dot(target_vec>=0, 1-score_vec)
        if den > 0:
            return num / den
        elif den == 0:
            return 0

    score_vec = np.zeros(matrix.shape[0])
    for offset in sorted(range(-matrix.DIST, matrix.DIST+1), key=lambda x: matrix.get_count_by_offset(x), reverse=True):
        target_vec = np.array(map(matrix.get_positions_by_offset(offset), lambda x: indicator(x, gene_name)))
        phat = caliculate_phat(target_vec, score_vec)
        target_vec[target_vec==-1] = phat #fill missing value with estimated phat
        score_vec += target_vec * (1 - score_vec)
    score = np.mean(score_vec)
    return score

