# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
import re
import warnings

from .variants import VariantRecord, BracketSVRecord, ShorthandSVRecord, VariantType

# Regex for SVs
BRACKET_SV_REGEX = re.compile(r'([.A-Za-z]*)(\[|\])([^\]\[:]+:[0-9]+)(\[|\])([.A-Za-z]*)')
SHORTHAND_SV_REGEX = re.compile(r'<(DEL|INS|DUP|INV|CNV])(:[A-Za-z]+)*>')
SGL_SV_REGEX = re.compile(r'\.[.A-Za-z]+|[.A-Za-z]+\.')
STANDARD_RECORD_REGEX = re.compile(r'([.A-Za-z]+)')


def _parse_bracket_sv(rec):
    sv_match_bracket = BRACKET_SV_REGEX.fullmatch(rec.alts[0])
    if not sv_match_bracket:
        return None
    # Extract ALT data from regex
    alt_prefix = sv_match_bracket.group(1)
    alt_bracket = sv_match_bracket.group(2)
    alt_contig, alt_pos = sv_match_bracket.group(3).split(':')
    alt_suffix = sv_match_bracket.group(5)
    alt_sv_bracket = BracketSVRecord(alt_prefix, alt_bracket, alt_contig, int(alt_pos), alt_suffix)
    # End position
    end_pos = int(alt_pos) if alt_contig == rec.contig else rec.stop

    # Extract type
    if rec.contig != alt_contig:
        # BND & INS with different contig ~ TRN
        variant_type = VariantType.TRN
        length = 0
    else:
        # INV -> 1 10 N]1:20] or 1 20 N]1:10]
        #        1 10 [1:20[N or 1 20 [1:10[N
        # DEL -> 1 10 N[1:20[ or 1 20 ]1:10]N
        # DUP -> 1 10 ]1:20]N or 1 20 N[1:10[
        length = abs(end_pos - rec.pos)
        equivalent_bracket = alt_bracket
        equivalent_preffix = alt_prefix
        # Transform REF/ALT to equivalent notation
        if int(alt_pos) < rec.pos:
            equivalent_bracket = ']' if alt_bracket == '[' else '['
            equivalent_preffix = alt_suffix
        if equivalent_preffix and equivalent_bracket == '[':
            variant_type = VariantType.DEL
        elif not equivalent_preffix and equivalent_bracket == ']':
            variant_type = VariantType.DUP
        else:
            variant_type = VariantType.INV

    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, end_pos, length, rec.id, rec.ref, rec.alts[0],
                               rec.qual, rec.filter, rec.info, variant_type, alt_sv_bracket, None)
    return vcf_record


def _parse_shorthand_sv(rec):
    sv_match_shorthand = SHORTHAND_SV_REGEX.fullmatch(rec.alts[0])
    if not sv_match_shorthand:
        return None
    # Extract ALT data from regex
    alt_type = sv_match_shorthand.group(1)
    alt_extra = sv_match_shorthand.group(2).split(':') if sv_match_shorthand.group(2) else None
    alt_sv_shorthand = ShorthandSVRecord(alt_type, alt_extra)
    length = abs(rec.stop - rec.pos)

    # Extract type
    if alt_type == 'DEL':
        variant_type = VariantType.DEL
    elif alt_type == 'INS':
        if not 'SVLEN' in rec.info:
            warnings.warn(f'SVLEN not found in INFO field for <INS> shorthand record. Defaults to 0.')
            length = 0
        elif isinstance(rec.info['SVLEN'], tuple):
            length = rec.info['SVLEN'][0]
        else:
            length = rec.info['SVLEN']
        length = abs(int(length))
        variant_type = VariantType.INS
    elif alt_type == 'DUP':
        variant_type = VariantType.DUP
    elif alt_type == 'INV':
        variant_type = VariantType.INV
    elif alt_type == 'CNV':
        variant_type = VariantType.CNV
    else:
         raise ValueError(f'Unknown variant type: {alt_type}. Skipping:\n{rec}')

    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, length, rec.id, rec.ref, rec.alts[0],
                               rec.qual, rec.filter, rec.info, variant_type, None, alt_sv_shorthand)
    return vcf_record


def _parse_sgl_sv(rec):
    sv_match_sgl = SGL_SV_REGEX.fullmatch(rec.alts[0])
    if not sv_match_sgl or 'SVTYPE' not in rec.info:
        return None
    variant_type = VariantType.SGL
    length = 0
    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, length, rec.id, rec.ref, rec.alts[0],
                               rec.qual, rec.filter, rec.info, variant_type, None, None)
    return vcf_record


def _parse_standard_record(rec):
    match = STANDARD_RECORD_REGEX.fullmatch(rec.alts[0])
    if not match:
        return None
    if len(rec.alts[0]) == len(rec.ref):
        length = 0
        variant_type = VariantType.SNV
    elif len(rec.alts[0]) > len(rec.ref):
        length = len(rec.alts[0]) - 1
        variant_type = VariantType.INS
    else:
        length = len(rec.ref) - 1
        variant_type = VariantType.DEL
    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, length, rec.id, rec.ref, rec.alts[0],
                               rec.qual, rec.filter, rec.info, variant_type, None, None)
    return vcf_record
