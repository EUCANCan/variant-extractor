# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
import re

from .variants import BracketSVRecord

NUMBER_CONTIG_REGEX = re.compile(r'[0-9]+')


def _select_record(variant_record_1, variant_record_2):
    # Same contig, select lowest position
    if variant_record_1.contig == variant_record_2.contig:
        return variant_record_1 if variant_record_1.pos < variant_record_2.pos else variant_record_2
    # Different contig
    else:
        match_1 = NUMBER_CONTIG_REGEX.search(variant_record_1.contig)
        match_2 = NUMBER_CONTIG_REGEX.search(variant_record_2.contig)
        # Both contigs do not contain numbers or follow different structure, select lowest in lexicographical order
        if not match_1 or not match_2 or match_1.start() != match_2.start():
            return variant_record_1 if variant_record_1.contig < variant_record_2.contig else variant_record_2
        else:
            # Both contigs contain numbers, select lowest number
            record = variant_record_1 if int(match_1.group()) < int(match_2.group()) else variant_record_2
            return record


def _permute_bracket_sv(variant_record):
    # Transform REF/ALT to equivalent notation so that REF contains the lowest position
    new_contig = variant_record.alt_sv_bracket.contig
    new_pos = variant_record.alt_sv_bracket.pos
    new_end = new_pos
    alt_prefix = variant_record.alt_sv_bracket.suffix
    alt_suffix = variant_record.alt_sv_bracket.prefix
    alt_bracket = ']' if variant_record.alt_sv_bracket.bracket == '[' else '['
    alt_contig = variant_record.contig
    alt_pos = variant_record.pos
    new_alt = f'{alt_prefix}{alt_bracket}{alt_contig}:{alt_pos}{alt_bracket}{alt_suffix}'
    alt_sv_bracket = BracketSVRecord(alt_prefix, alt_bracket, alt_contig, alt_pos, alt_suffix)
    variant_record = variant_record._replace(contig=new_contig, pos=new_pos,
                                             end=new_end, alt=new_alt, alt_sv_bracket=alt_sv_bracket)
    return variant_record


def _convert_inv_to_bracket(variant_record):
    # Convert INV to equivalent bracket notation. Ex:
    # 2 321682 T <INV> END=421681
    # is equivalent to:
    # 2 321681 . .]2:421681]
    # 2 321682 T [2:421682[T
    alt_1 = f'.]{variant_record.contig}:{variant_record.end}]'
    alt_sv_bracket_1 = BracketSVRecord('.', ']', variant_record.contig, variant_record.end, None)
    length_1 = abs(variant_record.end - (variant_record.pos - 1))
    variant_record_1 = variant_record._replace(
        pos=variant_record.pos-1, length=length_1, id=variant_record.id+'_1', ref='.', alt=alt_1, alt_sv_bracket=alt_sv_bracket_1)

    alt_2 = f'[{variant_record.contig}:{variant_record.end+1}[{variant_record.ref}'
    alt_sv_bracket_2 = BracketSVRecord(None, '[', variant_record.contig, variant_record.pos, variant_record.ref)
    length_2 = abs(variant_record.end + 1 - variant_record.pos)
    variant_record_2 = variant_record._replace(
        end=variant_record.end+1, length=length_2, id=variant_record.id+'_2', alt=alt_2, alt_sv_bracket=alt_sv_bracket_2)

    return variant_record_1, variant_record_2
