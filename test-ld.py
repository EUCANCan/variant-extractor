import numpy as np


def print_distance_matrix(a, b, matrix):
    matrix_string = '\t\t'
    matrix_string += '\t'.join(b) + '\n'
    for i in range(len(a) + 1):
        if i != 0:
            matrix_string += a[i-1]
        for j in range(len(b) + 1):
            matrix_string += f'\t{matrix[i][j]}'
        matrix_string += '\n'
    print(matrix_string)


def levenshtein_distance_matrix(a, b):
    matrix = np.zeros((len(a) + 1, len(b) + 1), dtype=int)
    matrix[:, 0] = np.arange(len(a) + 1)
    matrix[0] = np.arange(len(b) + 1)
    print_distance_matrix(a, b, matrix)

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i-1] == b[j-1]:
                matrix[i][j] = matrix[i-1][j-1]
            else:
                matrix[i][j] = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]) + 1
    print_distance_matrix(a, b, matrix)
    return matrix

class Variant(NamedTuple):
    chrom: str
    start_pos: int
    seq: str
    length: int
    start_base: str

def get_variants_from_matrix(ref, alt, matrix):
    end_pos = pos + len(ref)
    last_pos = end_pos
    variants = []
    i = len(ref)
    j = len(alt)
    while i > 0 or j > 0:
        current_cell = matrix[i][j]
        diagonal_cell = matrix[i-1][j-1] if i > 0 and j > 0 else current_cell + 1
        up_cell = matrix[i-1][j] if i > 0 else current_cell + 1
        left_cell = matrix[i][j-1] if j > 0 else current_cell + 1
        if diagonal_cell <= up_cell and diagonal_cell <= left_cell:
            if diagonal_cell < current_cell:
                ops.append(f'{i-1} Replace {ref[i-1]} with {b[j-1]}')
            i -= 1
            j -= 1
            last_pos = i
        elif left_cell <= up_cell and left_cell < current_cell:
            ops.append(f'{i-1} Insert {alt[j-1]}')
            j -= 1
        else:
            ops.append(f'{i-1} Delete {ref[i-1]}')
            i -= 1
    return ops[::-1]


if __name__ == '__main__':
    ref = 'GCT'
    alt = 'GGAGTAAGAAAC'
    distance_matrix = levenshtein_distance_matrix(ref, alt)
    ops = get_ops_from_matrix(ref, alt, distance_matrix)
    print(ops)
