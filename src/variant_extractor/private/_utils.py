# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# MIT License
import re

from ..variants import BreakendSVRecord, VariantRecord, VariantType

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


def permute_breakend_sv(variant_record: VariantRecord, fasta_ref=None):
    assert variant_record.alt_sv_breakend is not None
    # Transform REF/ALT to equivalent notation
    # Equivalencies:
    # 1 500 . N N[7:800[ 	7 800 . N ]1:500]N
    # 1 500 . N ]7:800]N 	7 800 . N N[1:500[
    # 1 500 . N [7:800[N 	7 800 . N [1:500[N
    # 1 500 . N N]7:800] 	7 800 . N N]1:500]
    new_contig = variant_record.alt_sv_breakend.contig
    alt_contig = variant_record.contig
    new_pos = variant_record.alt_sv_breakend.pos
    alt_pos = variant_record.pos
    new_end = alt_pos if new_contig == alt_contig else new_pos
    if variant_record.alt_sv_breakend.prefix and variant_record.alt_sv_breakend.bracket == '[':
        alt_prefix = None
        alt_suffix = variant_record.alt_sv_breakend.prefix if new_contig == alt_contig else 'N'
        if alt_suffix == 'N' and fasta_ref is not None:
            alt_suffix = fasta_ref.fetch(new_contig, new_pos-1, new_pos).upper()
        ref = alt_suffix
        alt_breakend = ']'
    elif variant_record.alt_sv_breakend.suffix and variant_record.alt_sv_breakend.bracket == ']':
        alt_prefix = variant_record.alt_sv_breakend.suffix if new_contig == alt_contig else 'N'
        if alt_prefix == 'N' and fasta_ref is not None:
            alt_prefix = fasta_ref.fetch(new_contig, new_pos-1, new_pos).upper()
        ref = alt_prefix
        alt_suffix = None
        alt_breakend = '['
    else:
        alt_prefix = variant_record.alt_sv_breakend.prefix
        alt_suffix = variant_record.alt_sv_breakend.suffix
        ref = variant_record.ref
        alt_breakend = variant_record.alt_sv_breakend.bracket
    new_alt = f'{alt_prefix if alt_prefix else ""}{alt_breakend}{alt_contig}:{alt_pos}{alt_breakend}{alt_suffix if alt_suffix else ""}'
    alt_sv_breakend = BreakendSVRecord(alt_prefix, alt_breakend, alt_contig, alt_pos, alt_suffix)
    variant_record = variant_record._replace(contig=new_contig, pos=new_pos, ref=ref,
                                             end=new_end, alt=new_alt, alt_sv_breakend=alt_sv_breakend)
    return variant_record


def convert_del_to_ins(variant_record: VariantRecord, fasta_ref=None):
    assert variant_record.alt_sv_breakend is not None
    # Convert DEL to INS
    # 1 100 T TATATATATACACAC[1:101[
    # 1 101 A ]1:100]ATATATATACACACA
    # 1 100 T TATATATATACACAC
    if variant_record.alt_sv_breakend.bracket == '[':
        assert variant_record.alt_sv_breakend.prefix is not None
        pos = variant_record.pos
        ref = variant_record.ref
        alt = variant_record.alt_sv_breakend.prefix
    else:
        assert variant_record.alt_sv_breakend.suffix is not None
        pos = variant_record.alt_sv_breakend.pos
        ref = 'N' if fasta_ref is None else fasta_ref.fetch(variant_record.contig, pos-1, pos).upper()
        alt = ref + variant_record.alt_sv_breakend.suffix
    length = len(alt) - 1
    return variant_record._replace(pos=pos, end=pos, ref=ref, alt=alt, length=length, alt_sv_breakend=None, variant_type=VariantType.INS)


def convert_inv_to_breakend(variant_record: VariantRecord, fasta_ref=None):
    # Convert INV to equivalent breakend notation. Ex:
    # 2 321682 T <INV> END=421681
    # is equivalent to:
    # 2 321681 . .]2:421681]
    # 2 321682 T [2:421682[T
    ref_1 = 'N' if fasta_ref is None else \
        fasta_ref.fetch(variant_record.contig, variant_record.pos-2, variant_record.pos - 1).upper()
    alt_1 = f'{ref_1}]{variant_record.contig}:{variant_record.end}]'
    alt_sv_breakend_1 = BreakendSVRecord(ref_1, ']', variant_record.contig, variant_record.end, None)
    length_1 = abs(variant_record.end - (variant_record.pos - 1))
    variant_record_1 = variant_record._replace(
        pos=variant_record.pos-1, length=length_1, id=variant_record.id+'_1' if variant_record.id else None, ref=ref_1, alt=alt_1, alt_sv_breakend=alt_sv_breakend_1, alt_sv_shorthand=None)

    alt_2 = f'[{variant_record.contig}:{variant_record.end+1}[{variant_record.ref}'
    alt_sv_breakend_2 = BreakendSVRecord(None, '[', variant_record.contig, variant_record.pos, variant_record.ref)
    length_2 = abs(variant_record.end + 1 - variant_record.pos)
    variant_record_2 = variant_record._replace(
        end=variant_record.end+1, length=length_2, id=variant_record.id+'_2' if variant_record.id else None, alt=alt_2, alt_sv_breakend=alt_sv_breakend_2, alt_sv_shorthand=None)

    return variant_record_1, variant_record_2
