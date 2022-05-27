import numpy as np
from typing import NamedTuple
from enum import Enum

class Inertia(Enum):
    DIAGONAL = 1
    RIGHT = 2
    BOTTOM = 3

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

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i-1] == b[j-1]:
                matrix[i][j] = matrix[i-1][j-1]
            else:
                matrix[i][j] = min(matrix[i-1][j-1], matrix[i-1][j], matrix[i][j-1]) + 1
    return matrix

class Variant(NamedTuple):
    chrom: str
    start_pos: int
    seq: str
    length: int
    start_base: str

def handle_ins(chrom, pos, start_base, seq):
    return Variant(chrom, pos, seq, len(seq), start_base)

def handle_del(chrom, pos, start_base, seq):
    return Variant(chrom, pos, seq, -len(seq), start_base)

def handle_snv(chrom, pos, base_ref, base_alt):
    return Variant(chrom, pos, base_alt, 0, base_ref)

def get_variants_from_matrix(ref, alt, matrix, chrom, pos):
    variants = []
    var_type = ''
    i = 0
    j = 0
    start_pos_ref = i
    start_pos_alt = j

    def prev_base(i):
        return ref[i+1] if i < len(ref) else 'N'

    inertia = Inertia.DIAGONAL
    while i < len(ref) or j < len(alt):
        current_cell = matrix[i][j]
        diagonal_cell = matrix[i+1][j+1] if i < len(ref) and j < len(alt) else current_cell + 99
        bottom_cell = matrix[i+1][j] if i < len(ref) else current_cell + 99
        right_cell = matrix[i][j+1] if j < len(alt) else current_cell + 99

        diagonal_is_walkable = diagonal_cell <= bottom_cell and diagonal_cell <= right_cell
        bottom_is_walkable = bottom_cell <= diagonal_cell and bottom_cell <= right_cell
        right_is_walkable = right_cell <= diagonal_cell and right_cell <= bottom_cell

        print(f'({i},{j}) {inertia} {current_cell}')

        if (inertia == Inertia.DIAGONAL and diagonal_is_walkable) or (not bottom_is_walkable and not right_is_walkable):
            if var_type == '+':
                variants.append(handle_ins(chrom, pos+i-1, ref[i-1], alt[start_pos_alt:j]))
            elif var_type == '-':
                variants.append(handle_del(chrom, pos+start_pos_ref-1, ref[start_pos_ref-1], ref[start_pos_ref:i]))
            if diagonal_cell != current_cell:
                variants.append(handle_snv(chrom, pos+i, ref[i], alt[j]))
            i += 1
            j += 1
            var_type = ''
            inertia = Inertia.DIAGONAL
        elif (inertia == Inertia.RIGHT and right_is_walkable) or not bottom_is_walkable:
            # Insertion
            if var_type == '-':
                variants.append(handle_del(chrom, pos+start_pos_ref-1, ref[start_pos_ref-1], ref[start_pos_ref:i]))
            if var_type != '+':
                start_pos_alt = j
            var_type = '+'
            j += 1
            inertia = Inertia.RIGHT
        else:
            # Deletion
            if var_type == '+':
                variants.append(handle_ins(chrom, pos+i-1, ref[i-1], alt[start_pos_alt:j]))
            if var_type != '-':
                start_pos_ref = i
            var_type = '-'
            i += 1
            inertia = Inertia.BOTTOM
    if var_type == '+':
        variants.append(handle_ins(chrom, pos+i-1, ref[i-1], alt[start_pos_alt:j]))
    elif var_type == '-':
        variants.append(handle_del(chrom, pos+start_pos_ref-1, ref[start_pos_ref-1], ref[start_pos_ref:i]))
    return variants



def variant_str(var):
    if var.length == 0:
        return f'{var.chrom}\t{var.start_pos}\t{var.start_base}\t{var.seq}'
    elif var.length > 0:
        return f'{var.chrom}\t{var.start_pos}\t{var.start_base}\t{var.start_base}{var.seq}'
    else:
        return f'{var.chrom}\t{var.start_pos}\t{var.start_base}{var.seq}\t{var.start_base}'

if __name__ == '__main__':
    ref = 'AAGAGTGCGCATAAGATGAGGATATGG'
    alt = 'ACATTTAT'
    ref = ref.upper()
    alt = alt.upper()
    chrom = '10'
    pos = 111967168
    distance_matrix = levenshtein_distance_matrix(ref, alt)
    print_distance_matrix(ref, alt, distance_matrix)
    variants = get_variants_from_matrix(ref, alt, distance_matrix, chrom, pos)
    for var in variants:
        print(variant_str(var))
