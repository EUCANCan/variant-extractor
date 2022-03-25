# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# MIT License
import re

from .variants import VariantRecord, BracketSVRecord, ShorthandSVRecord

# Regex for SVs
BRACKET_SV_REGEX = re.compile(r'([.A-Za-z]*)(\[|\])([^\]\[:]+:[0-9]+)(\[|\])([.A-Za-z]*)')
SHORTHAND_SV_REGEX = re.compile(r'<(DEL|INS|DUP|INV|CNV])(:[A-Za-z]+)*>')
SGL_SV_REGEX = re.compile(r'\.?[.A-Za-z]\.?')

NUMBER_CONTIG_REGEX = re.compile(r'[0-9]+')


def select_record(variant_record_1, variant_record_2):
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


def permute_bracket_sv(variant_record):
    # Transform REF/ALT to equivalent notation so that REF contains the lowest position
    new_contig = variant_record.alt_sv_bracket.contig
    new_pos = variant_record.alt_sv_bracket.pos
    new_end = new_pos
    alt_prefix = variant_record.alt_sv_bracket.suffix
    alt_suffix = variant_record.alt_sv_bracket.prefix
    alt_bracket = ']' if variant_record.alt_sv_bracket.bracket == '[' else '['
    alt_contig = variant_record.contig
    alt_pos = variant_record.pos
    new_alts = [f'{alt_prefix}{alt_bracket}{alt_contig}:{alt_pos}{alt_bracket}{alt_suffix}']
    alt_sv_bracket = BracketSVRecord(alt_prefix, alt_bracket, alt_contig, alt_pos, alt_suffix)
    variant_record = VariantRecord(new_contig, new_pos, new_end, variant_record.id, variant_record.ref, new_alts,
                                variant_record.filter, variant_record.info, alt_sv_bracket, None)
    return variant_record


def extract_bracket_sv(rec):
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


def extract_shorthand_sv(rec):
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


def extract_sgl_sv(rec):
    sv_match_sgl = SGL_SV_REGEX.search(rec.alts[0])
    if not sv_match_sgl or 'SVTYPE' not in rec.info:
        return None
    # Create new record
    vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, rec.id, rec.ref,
                               rec.alts, rec.filter, rec.info, None, None)
    return vcf_record
