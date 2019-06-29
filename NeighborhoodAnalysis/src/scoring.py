import numpy as np

def score_naive(gene_name, matrix):
    total = matrix.shape[0]
    found = len(set([pos.i for pos in matrix.get_positions_by_gene_name(gene_name)]))
    if total > 0:
        return found / total
    else:
        return 0

def score_independent(gene_name, matrix):
    prod = 1
    for offset in range(-matrix.DIST, matrix.DIST+1):
        total = matrix.get_count_by_offset(offset)
        if total > 0:
            found = sum([pos.gene_name==gene_name for pos in matrix.get_positions_by_offset(offset)])
            prod *= 1 - found / total
    return 1 - prod

def score_conditional(gene_name, matrix):

	def indicator(position, gene_name):
		if position is None:
			return -1.0
		elif position.gene_name != gene_name:
			return 0.0
		else:
			return 1.0

    def caliculate_phat(target_vec, score_vec):
        num = np.dot(target_vec==1, 1-score_vec)
        den = np.dot(target_vec>=0, 1-score_vec)
        assert num >= 0 and den >= 0
        if den > 0:
            return num / den
        elif den == 0:
            return 0

	score_vec = np.zeros(matrix.shape[0])
	offsets = list(range(-matrix.DIST, 0)) + list(range(1, matrix.DIST+1))
	for offset in sorted(offsets, key=lambda x: matrix.get_count_by_offset(x), reverse=True):
		target_vec = np.array([indicator(pos, gene_name) for pos in matrix.get_positions_by_offset(offset)])
		phat = caliculate_phat(target_vec, score_vec)
		target_vec[target_vec==-1] = phat #fill missing value with estimated phat
		score_vec += target_vec * (1 - score_vec)
    score = np.mean(score_vec)
    return score

