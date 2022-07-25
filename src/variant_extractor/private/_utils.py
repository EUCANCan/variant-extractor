# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
import re

from ..variants import BracketSVRecord, VariantRecord, VariantType

NUMBER_CONTIG_REGEX = re.compile(r'[0-9]+')


def compare_contigs(contig_1, contig_2):
    # Same contig
    if contig_1 == contig_2:
        return 0
    match_1 = NUMBER_CONTIG_REGEX.search(contig_1)
    match_2 = NUMBER_CONTIG_REGEX.search(contig_2)
    # Both contigs do not contain numbers or follow different structure, select lowest in lexicographical order
    if not match_1 or not match_2 or match_1.start() != match_2.start():
        return -1 if contig_1 <= contig_2 else 1
    else:
        # Both contigs contain numbers, select lowest number
        return -1 if int(match_1.group()) <= int(match_2.group()) else 1


def permute_bracket_sv(variant_record: VariantRecord, fasta_ref=None):
    # Transform REF/ALT to equivalent notation
    # Equivalencies:
    # 1 500 . N N[7:800[ 	7 800 . N ]1:500]N
    # 1 500 . N ]7:800]N 	7 800 . N N[1:500[
    # 1 500 . N [7:800[N 	7 800 . N [1:500[N
    # 1 500 . N N]7:800] 	7 800 . N N]1:500]
    new_contig = variant_record.alt_sv_bracket.contig
    alt_contig = variant_record.contig
    new_pos = variant_record.alt_sv_bracket.pos
    alt_pos = variant_record.pos
    new_end = alt_pos if new_contig == alt_contig else new_pos
    if variant_record.alt_sv_bracket.prefix and variant_record.alt_sv_bracket.bracket == '[':
        alt_prefix = None
        alt_suffix = variant_record.alt_sv_bracket.prefix if new_contig == alt_contig else 'N'
        if alt_suffix == 'N' and fasta_ref is not None:
            alt_suffix = fasta_ref.fetch(new_contig, new_pos-1, new_pos).upper()
        ref = alt_suffix
        alt_bracket = ']'
    elif variant_record.alt_sv_bracket.suffix and variant_record.alt_sv_bracket.bracket == ']':
        alt_prefix = variant_record.alt_sv_bracket.suffix if new_contig == alt_contig else 'N'
        if alt_prefix == 'N' and fasta_ref is not None:
            alt_prefix = fasta_ref.fetch(new_contig, new_pos-1, new_pos).upper()
        ref = alt_prefix
        alt_suffix = None
        alt_bracket = '['
    else:
        alt_prefix = variant_record.alt_sv_bracket.prefix
        alt_suffix = variant_record.alt_sv_bracket.suffix
        ref = variant_record.ref
        alt_bracket = variant_record.alt_sv_bracket.bracket
    new_alt = f'{alt_prefix if alt_prefix else ""}{alt_bracket}{alt_contig}:{alt_pos}{alt_bracket}{alt_suffix if alt_suffix else ""}'
    alt_sv_bracket = BracketSVRecord(alt_prefix, alt_bracket, alt_contig, alt_pos, alt_suffix)
    variant_record = variant_record._replace(contig=new_contig, pos=new_pos, ref=ref,
                                             end=new_end, alt=new_alt, alt_sv_bracket=alt_sv_bracket)
    return variant_record


def convert_del_to_ins(variant_record: VariantRecord, fasta_ref=None):
    # Convert DEL to INS
    # 1 100 T TATATATATACACAC[1:101[
    # 1 101 A ]1:100]ATATATATACACACA
    # 1 100 T TATATATATACACAC
    if variant_record.alt_sv_bracket.bracket == '[':
        pos = variant_record.pos
        ref = variant_record.ref
        alt = variant_record.alt_sv_bracket.prefix
    else:
        pos = variant_record.alt_sv_bracket.pos
        ref = 'N' if fasta_ref is None else fasta_ref.fetch(variant_record.contig, pos-1, pos).upper()
        alt = ref + variant_record.alt_sv_bracket.suffix
    length = len(alt) - 1
    return variant_record._replace(pos=pos, end=pos, ref=ref, alt=alt, length=length, alt_sv_bracket=None, variant_type=VariantType.INS)


def convert_inv_to_bracket(variant_record: VariantRecord, fasta_ref=None):
    # Convert INV to equivalent bracket notation. Ex:
    # 2 321682 T <INV> END=421681
    # is equivalent to:
    # 2 321681 . .]2:421681]
    # 2 321682 T [2:421682[T
    ref_1 = 'N' if fasta_ref is None else \
        fasta_ref.fetch(variant_record.contig, variant_record.pos-2, variant_record.pos - 1).upper()
    alt_1 = f'{ref_1}]{variant_record.contig}:{variant_record.end}]'
    alt_sv_bracket_1 = BracketSVRecord(ref_1, ']', variant_record.contig, variant_record.end, None)
    length_1 = abs(variant_record.end - (variant_record.pos - 1))
    variant_record_1 = variant_record._replace(
        pos=variant_record.pos-1, length=length_1, id=variant_record.id+'_1' if variant_record.id else None, ref=ref_1, alt=alt_1, alt_sv_bracket=alt_sv_bracket_1, alt_sv_shorthand=None)

    alt_2 = f'[{variant_record.contig}:{variant_record.end+1}[{variant_record.ref}'
    alt_sv_bracket_2 = BracketSVRecord(None, '[', variant_record.contig, variant_record.pos, variant_record.ref)
    length_2 = abs(variant_record.end + 1 - variant_record.pos)
    variant_record_2 = variant_record._replace(
        end=variant_record.end+1, length=length_2, id=variant_record.id+'_2' if variant_record.id else None, alt=alt_2, alt_sv_bracket=alt_sv_bracket_2, alt_sv_shorthand=None)

    return variant_record_1, variant_record_2
