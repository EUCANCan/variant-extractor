import pysam
import difflib
from typing import NamedTuple


class Variant(NamedTuple):
    chrom: str
    start_pos: int
    seq: str
    length: int
    start_base: str


def find_differences(chrom, pos, padding, old, new):
    print(old)
    print(new)
    diff = difflib.ndiff(old.upper(), new.upper())
    print('|'.join(diff))
    diff = difflib.ndiff(old.upper(), new.upper())

    variant_pos = pos - padding - 1
    building_seq = ''
    found_variants = []
    seq_type = ' '
    last_char = ''

    for line in diff:
        start_char, _, char = line
        if start_char == ' ':
            if len(building_seq) > 0:
                var_len = len(building_seq) if seq_type == '+' else -len(building_seq)
                found_variants.append(Variant(chrom, variant_pos, building_seq, var_len, last_char))
                if seq_type == '-':
                    variant_pos += len(building_seq)
                building_seq = ''
                seq_type = ' '
            last_char = char
            variant_pos += 1
        elif start_char == '-':
            if seq_type == '+':
                found_variants.append(Variant(chrom, variant_pos, building_seq, len(building_seq), last_char))
                building_seq = ''
            building_seq += char
            seq_type = start_char
        elif start_char == '+':
            if seq_type == '-':
                found_variants.append(Variant(chrom, variant_pos, building_seq, -len(building_seq), last_char))
                building_seq = ''
                variant_pos += len(building_seq)
            building_seq += char
            seq_type = start_char
        else:
            raise Exception(f'Unknown diff line: {line}')

    if len(building_seq) > 0:
        var_len = len(building_seq) if seq_type == '+' else -len(building_seq)
        found_variants.append(Variant(chrom, variant_pos, building_seq, var_len, last_char))

    return found_variants


class VariantNormalizator:
    def __init__(self, fasta_ref):
        self.fasta = pysam.FastaFile(fasta_ref)

    def normalize(self, chrom, pos, ref, alt):
        padding = 50
        ref = ref.upper()
        alt = alt.upper()
        before = self.fasta.fetch(chrom, pos-padding-1, pos-1).lower()
        middle = self.fasta.fetch(chrom, pos-1, pos+len(ref)-1).upper()
        after = self.fasta.fetch(chrom, pos+len(ref)-1, pos+len(ref)+padding-1).lower()
        assert middle == ref, f'Reference mismatch: {middle} != {ref}'
        old = before + middle + after
        new = before + alt + after
        found_variants = find_differences(chrom, pos, padding, old, new)
        variants_str = self.build_variants(chrom, pos, padding, old, new, found_variants)
        return variants_str

    def build_variants(self, chrom, pos, padding, old, new, found_variants):
        # Search for SNVs
        i = 0
        req_offset = 0
        alt_offset = 0
        str_offset = pos - padding
        # print(f'str offset {str_offset}')
        # print('old len =', len(old))
        # print('new len =', len(new))
        variants_str = []
        non_matched_variants = []
        while i < len(found_variants):
            matched = False
            curr_variant = found_variants[i]
            size_diff = curr_variant.length
            for j in range(i+1, len(found_variants)):
                next_variant = found_variants[j]
                size_diff += next_variant.length
                if size_diff == 0:
                    # print(f'SNVS {found_variants[i:j+1]}')
                    var_start_pos = curr_variant.start_pos + 1
                    var_end_pos = next_variant.start_pos+1
                    if next_variant.length < 0:
                        var_end_pos -= next_variant.length
                    elif curr_variant.start_pos == next_variant.start_pos and curr_variant.length < 0:
                        var_end_pos -= curr_variant.length
                    ref_str_start = var_start_pos-str_offset
                    ref_str_end = var_end_pos-str_offset
                    alt_str_start = alt_offset+var_start_pos-str_offset
                    alt_str_end = alt_offset+var_end_pos-str_offset
                    ref_seq = old[ref_str_start:ref_str_end]
                    alt_seq = new[alt_str_start:alt_str_end]
                    # print(ref_str_start, ref_str_end, ref_seq, var_start_pos, var_end_pos)
                    # print(alt_str_start, alt_str_end, alt_seq, var_start_pos, var_end_pos)
                    variants_str += self.handle_snv(chrom, var_start_pos, ref_seq, alt_seq)
                    i = j + 1
                    matched = True
                    break
            if not matched:
                non_matched_variants.append(curr_variant)
                alt_offset += curr_variant.length
                i += 1
        variants_str += self.handle_non_matched_variants(non_matched_variants)
        return variants_str

    def handle_snv(self, chrom, start_pos, ref_seq, alt_seq):
        # print('SNV')
        variants_str = []
        for pos in range(len(ref_seq)):
            ref_base = ref_seq[pos]
            alt_base = alt_seq[pos]
            if ref_base != alt_base:
                variants_str.append(f'{chrom}\t{start_pos}\t{ref_base}\t{alt_base}')
            start_pos += 1
        return variants_str

    def handle_non_matched_variants(self, non_matched_variants):
        # print('REST')
        variants_str = []
        i = 0
        while i < len(non_matched_variants):
            if i < len(non_matched_variants) - 1 and non_matched_variants[i].start_pos == non_matched_variants[i+1].start_pos:
                variants_str += self.handle_compound_indel(non_matched_variants[i], non_matched_variants[i+1])
                i += 2
            else:
                variants_str.append(self.handle_standard_variant(non_matched_variants[i]))
                i += 1
        return variants_str

    def handle_standard_variant(self, variant):
        # Left trim
        while variant.start_base == variant.seq[-1]:
            new_pos = variant.start_pos-1
            new_start_base = self.fasta.fetch(variant.chrom, new_pos-1, new_pos).upper()
            variant = variant._replace(start_pos=variant.start_pos-1, seq=variant.start_base +
                                       variant.seq[:-1], start_base=new_start_base)

        chrom = variant.chrom
        pos = variant.start_pos
        if variant.length > 0:
            ref = variant.start_base
            alt = ref + variant.seq
        else:
            ref = variant.start_base + variant.seq
            alt = variant.start_base
        return f'{chrom}\t{pos}\t{ref}\t{alt}'

    def handle_compound_indel(self, variant_1, variant_2):
        # TODO: Trim?
        variants_str = []
        ref_first = variant_1.start_base
        del_variant, ins_variant = (variant_1, variant_2) if variant_1.length < 0 else (variant_2, variant_1)
        min_len = min(abs(del_variant.length), abs(ins_variant.length))
        pos = 0
        while pos < min_len:
            ref_base = del_variant.seq[pos]
            alt_base = ins_variant.seq[pos]
            if ref_base != alt_base:
                variants_str.append(f'{del_variant.chrom}\t{del_variant.start_pos+pos+1}\t{ref_base}\t{alt_base}')
            pos += 1
        if abs(del_variant.length) > abs(ins_variant.length):
            # Deletion
            variants_str.append(
                f'{del_variant.chrom}\t{del_variant.start_pos+pos}\t{del_variant.seq[pos-1:]}\t{del_variant.seq[pos-1]}')
        else:
            # Insertion
            variants_str.append(
                f'{ins_variant.chrom}\t{ins_variant.start_pos+pos}\t{del_variant.seq[pos-1]}\t{ins_variant.seq[pos-1:]}')
        return variants_str


if __name__ == '__main__':
    REF_FASTA = '../insilico-builder/test_data/ref_38/hg38.fa'
    VCF_INPUT = 'tests/only_pass.vcf'
    RESULT_VCF = 'tests/result.vcf'

    normalizator = VariantNormalizator(REF_FASTA)

    f_in = open(VCF_INPUT, 'r')
    f_out = open(RESULT_VCF, 'w')

    for line in f_in:
        chrom, pos, _, ref, alt = line.strip().split('\t')[:5]
        chrom = chrom.replace('chr', '')
        pos = int(pos)
        variants = normalizator.normalize(chrom, pos, ref, alt)
        for variant in variants:
            f_out.write(f'{variant}\n')
