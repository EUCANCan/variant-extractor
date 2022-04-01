# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
import re

from .variants import VariantRecord, BracketSVRecord, ShorthandSVRecord

# Regex for SVs
BRACKET_SV_REGEX = re.compile(r'([.A-Za-z]*)(\[|\])([^\]\[:]+:[0-9]+)(\[|\])([.A-Za-z]*)')
SHORTHAND_SV_REGEX = re.compile(r'<(DEL|INS|DUP|INV|CNV])(:[A-Za-z]+)*>')
SGL_SV_REGEX = re.compile(r'\.[.A-Za-z]+|[.A-Za-z]+\.')
STANDARD_RECORD_REGEX = re.compile(r'([.A-Za-z]+)')


def _parse_bracket_sv(rec):
    sv_match_bracket = BRACKET_SV_REGEX.search(rec.alts[0])
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
    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, end_pos, rec.id, rec.ref,
                               rec.alts, rec.filter, rec.info, alt_sv_bracket, None)
    return vcf_record


def _parse_shorthand_sv(rec):
    sv_match_shorthand = SHORTHAND_SV_REGEX.search(rec.alts[0])
    if not sv_match_shorthand:
        return None
    # Extract ALT data from regex
    alt_type = sv_match_shorthand.group(1)
    alt_extra = sv_match_shorthand.group(2).split(':') if sv_match_shorthand.group(2) else None
    alt_sv_shorthand = ShorthandSVRecord(alt_type, alt_extra)

    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, rec.id, rec.ref,
                               rec.alts, rec.filter, rec.info, None, alt_sv_shorthand)
    return vcf_record


def _parse_sgl_sv(rec):
    sv_match_sgl = SGL_SV_REGEX.search(rec.alts[0])
    if not sv_match_sgl or 'SVTYPE' not in rec.info:
        return None
    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, rec.id, rec.ref,
                               rec.alts, rec.filter, rec.info, None, None)
    return vcf_record


def _parse_standard_record(rec):
    match = STANDARD_RECORD_REGEX.search(rec.alts[0])
    if not match:
        return None
    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, rec.id, rec.ref,
                               rec.alts, rec.filter, rec.info, None, None)
    return vcf_record
