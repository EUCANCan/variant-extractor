import numpy as np
from typing import NamedTuple


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

def handle_ins(chrom, pos, start_base, seq):
    return Variant(chrom, pos, seq, len(seq), start_base)

def handle_del(chrom, pos, start_base, seq):
    return Variant(chrom, pos, seq, -len(seq), start_base)

def handle_snv(chrom, pos, base_ref, base_alt):
    return Variant(chrom, pos, base_alt, 0, base_ref)

def get_variants_from_matrix(ref, alt, matrix, chrom, pos):
    variants = []
    var_type = ''
    i = len(ref)
    j = len(alt)
    end_pos_ref = i
    end_pos_alt = j

    def prev_base(i):
        return ref[i-1] if i > 0 else 'N'

    while i > 0 or j > 0:
        current_cell = matrix[i][j]
        diagonal_cell = matrix[i-1][j-1] if i > 0 and j > 0 else current_cell + 1
        up_cell = matrix[i-1][j] if i > 0 else current_cell + 1
        left_cell = matrix[i][j-1] if j > 0 else current_cell + 1
        if diagonal_cell <= up_cell and diagonal_cell <= left_cell:
            if var_type == '+':
                variants.append(handle_ins(chrom, pos+i-1, prev_base(i), alt[j:end_pos_alt]))
            elif var_type == '-':
                variants.append(handle_del(chrom, pos+i-1, prev_base(i), ref[i:end_pos_ref]))
            if diagonal_cell < current_cell:
                variants.append(handle_snv(chrom, pos+i-1, prev_base(i), alt[j-1]))
            i -= 1
            j -= 1
            var_type = ''
        elif left_cell <= up_cell and left_cell < current_cell:
            # Insertion
            if var_type == '-':
                variants.append(handle_del(chrom, pos+i-1, prev_base(i), ref[i:end_pos_ref]))
            elif var_type != '+':
                end_pos_alt = j
            var_type = '+'
            j -= 1
        else:
            # Deletion
            if var_type == '+':
                variants.append(handle_ins(chrom, pos+i-1, prev_base(i), alt[j:end_pos_alt]))
            elif var_type != '-':
                end_pos_ref = i
            var_type = '-'
            i -= 1
    if var_type == '+':
        variants.append(handle_ins(chrom, pos+i-1, prev_base(i), alt[j:end_pos_alt]))
    elif var_type == '-':
        variants.append(handle_del(chrom, pos+i-1, prev_base(i), ref[i:end_pos_ref]))
    return variants[::-1]



def variant_str(var):
    if var.length == 0:
        return f'{var.chrom}\t{var.start_pos}\t{var.start_base}\t{var.seq}'
    elif var.length > 0:
        return f'{var.chrom}\t{var.start_pos}\t{var.start_base}\t{var.start_base}{var.seq}'
    else:
        return f'{var.chrom}\t{var.start_pos}\t{var.start_base}{var.seq}\t{var.start_base}'

if __name__ == '__main__':
    ref = 'GATGGAGACATG'
    alt = 'GTGT'
    ref = ref.upper()
    alt = alt.upper()
    chrom = '15'
    pos = 62309661
    distance_matrix = levenshtein_distance_matrix(ref, alt)
    variants = get_variants_from_matrix(ref, alt, distance_matrix, chrom, pos)
    for var in variants:
        print(variant_str(var))
