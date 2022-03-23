# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# MIT License
from typing import NamedTuple, Optional, List
from enum import Enum, auto
from os import path
from argparse import ArgumentParser
import math
import warnings
import re
import pysam

# TODO: What to do with REF=N ALT=.N (Single breakends)
# Regex for SVs
BRACKET_SV_REGEX = re.compile(r'([.A-Za-z]*)(\[|\])([^\]\[:]+:[0-9]+)(\[|\])([.A-Za-z]*)')
SHORTHAND_SV_REGEX = re.compile(r'<(DEL|INS|DUP|INV|CNV])(:[A-Z]+)*>')
NUMBER_CONTIG_REGEX = re.compile(r'[0-9]+')


class VariantType(Enum):
    """Enumeration with the different types of variations
    """
    SNV = auto()
    INDEL_INS = auto()
    INDEL_DEL = auto()
    DEL = auto()
    INS = auto()
    DUP = auto()
    INV = auto()
    CNV = auto()
    TRN = auto()


class BracketSVRecord(NamedTuple):
    """NamedTuple with the information of a bracketed SV record
    """
    prefix: Optional[str]
    """Prefix of the SV record with bracket notation. For example, for :code:`G]17:198982]` the prefix will be :code:`G`"""
    bracket: str
    """Bracket of the SV record with bracket notation. For example, for :code:`G]17:198982]` the bracket will be :code:`]`"""
    contig: str
    """Contig of the SV record with bracket notation. For example, for :code:`G]17:198982]` the contig will be :code:`17`"""
    pos: int
    """Position of the SV record with bracket notation. For example, for :code:`G]17:198982]` the position will be :code:`198982`"""
    suffix: Optional[str]
    """Suffix of the SV record with bracket notation. For example, for :code:`G]17:198982]` the suffix will be :code:`None`"""


class ShorthandSVRecord(NamedTuple):
    """NamedTuple with the information of a shorthand SV record
    """
    type: str
    """One of the following, :code:`'DEL'`, :code:`'INS'`, :code:`'DUP'`, :code:`'INV'` or :code:`'CNV'`"""
    extra: List[str]
    """Extra information of the SV. For example, for :code:`<DUP:TANDEM:AA>` the extra will be :code:`['TANDEM', 'AA']`"""


class VariantRecord(NamedTuple):
    """NamedTuple with the information of a variant record
    """
    contig: str
    """Contig name"""
    pos: int
    """Position of the variant in the contig"""
    end: int
    """End position of the variant in the contig (same as `pos` for TRN and SNV)"""
    id: str
    """Record identifier"""
    ref: str
    """Reference sequence"""
    alts: List[str]
    """List of alternative sequences"""
    filter: str
    """Original record filter"""
    info: dict
    """Original record info"""
    alt_sv_bracket: Optional[BracketSVRecord]
    """Bracketed SV info, present only for SVs with bracket notation. For example, :code:`G]17:198982]`"""
    alt_sv_shorthand: Optional[ShorthandSVRecord]
    """Shorthand SV info, present only for SVs with shorthand notation. For example, :code:`<DUP:TANDEM>`"""


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


def _extract_sv_from_brackets(rec):
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


def _extract_sv_from_shorthand(rec):
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


class VariantExtractor:
    """
    Reads and extracts variants from VCF files. This class is designed to be
    used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis.
    """

    def __init__(self, ensure_pairs=True):
        """
        Parameters
        ----------
        ensure_pairs : bool, optional
            If `True`, throws an exception if a breakend is missing a pair when all other were paired successfully.
        """
        self.ensure_pairs = ensure_pairs

    def read_vcf(self, vcf_file):
        """Reads VCF file and extracts all variants.

        Parameters
        ----------
        vcf_file : str
            A VCF formatted file. The file is automatically opened.

        Returns
        -------
        list[tuple[VariantType, VariantRecord]
        """

        self.__variants = []
        self.__pairs_found = 0
        self.__pending_sv_pairs = {}
        # Read the file
        with open(file=vcf_file, mode='r') as vcf_handle:
            save = pysam.set_verbosity(0)
            with pysam.VariantFile(vcf_handle) as vcf:
                pysam.set_verbosity(save)
                for rec in vcf:
                    self.__parse_record(rec)
        # Handle uncertain unpaired SV
        self.__handle_uncertain_sv()
        # Only single-paired records or not ensuring pairs
        if not self.ensure_pairs or self.__pairs_found == 0:
            if self.__pairs_found == 0:
                warnings.warn('WARNING: No SV breakend pairs found, assuming all records are single-paired')
            for alt_dicts in self.__pending_sv_pairs.values():
                for vcf_record in alt_dicts.values():
                    self.__parse_bracket_individual_sv(vcf_record)
        # Found unpaired records
        elif len(self.__pending_sv_pairs) > 0:
            exception_text = ''
            for alt_dicts in self.__pending_sv_pairs.values():
                for vcf_record in alt_dicts.values():
                    exception_text += f'REF={vcf_record.contig}:{vcf_record.pos}\tALT={vcf_record.alt_sv_bracket.contig}:{vcf_record.alt_sv_bracket.pos}\n'
            raise Exception(
                (f'Unpaired SV breakends:\n{exception_text}'
                 f'Exception: There are {len(self.__pending_sv_pairs)} unpaired SV breakends. '
                 'Please, check the entires show above in VCF file. '
                 'Use ensure_pairs=False to ignore unpaired SV breakends.'))
        return self.__variants

    def __parse_record(self, rec):
        if len(rec.alts) != 1:
            warnings.warn(f'WARNING: Skipping record with multiple alternate alleles ({rec})')
            return
        vcf_record = _extract_sv_from_brackets(rec)
        # Check if bracket SV record
        if vcf_record:
            return self.__parse_bracket_sv(vcf_record)
        # Check if shorthand SV record
        vcf_record = _extract_sv_from_shorthand(rec)
        if vcf_record:
            return self.__parse_shorthand_sv(vcf_record)
        # Indel or SNV
        vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, rec.id, rec.ref,
                                   rec.alts, rec.filter, rec.info, None, None)
        # Check if SNV
        if len(vcf_record.alts[0]) == len(vcf_record.ref):
            # REF=CTT ALT=ATG -> Normalize to 3 SNVs
            for i in range(len(vcf_record.ref)):
                if vcf_record.alts[0][i] != vcf_record.ref[i]:
                    new_vcf_record = vcf_record._replace(
                        ref=vcf_record.ref[i], pos=i+vcf_record.pos, end=i+vcf_record.end, alts=[vcf_record.alts[0][i]])
                    self.__variants.append((VariantType.SNV, new_vcf_record))
        # Check if INDEL_DEL
        elif len(vcf_record.alts[0]) < len(vcf_record.ref):
            return self.__variants.append((VariantType.INDEL_DEL, vcf_record))
        # Check if INDEL_INS
        else:
            return self.__variants.append((VariantType.INDEL_INS, vcf_record))

    def __parse_bracket_sv(self, vcf_record):
        # Check for pending SVs
        previous_record = self.__pop_pending_sv_pair(vcf_record)
        if previous_record is None:
            self.__store_pending_sv_pair(vcf_record)
            return
        # Mate SV found, parse it
        self.__pairs_found += 1
        record = _select_record(previous_record, vcf_record)
        self.__parse_bracket_individual_sv(record)

    def __store_pending_sv_pair(self, vcf_record):
        alt_name = f'{vcf_record.alt_sv_bracket.contig}{vcf_record.alt_sv_bracket.pos}'
        sv_name = f'{vcf_record.contig}{vcf_record.pos}'
        # Check if alt is already in the dictionary
        previous_records = self.__pending_sv_pairs.get(alt_name)
        if previous_records is None:
            # Create new entry for alt
            new_records = {sv_name: vcf_record}
            self.__pending_sv_pairs[alt_name] = new_records
        else:
            # Already exists alt entry, add new entry to alt
            previous_records[sv_name] = vcf_record

    def __pop_pending_sv_pair(self, vcf_record):
        sv_name = f'{vcf_record.contig}{vcf_record.pos}'
        previous_records = self.__pending_sv_pairs.get(sv_name)
        if previous_records is None:
            return None
        previous_record_alt_name = f'{vcf_record.alt_sv_bracket.contig}{vcf_record.alt_sv_bracket.pos}'
        previous_record_alt = previous_records.get(previous_record_alt_name)
        if previous_record_alt is None:
            return None
        previous_records.pop(previous_record_alt_name)
        if len(previous_records) == 0:
            self.__pending_sv_pairs.pop(sv_name)
        return previous_record_alt

    def __parse_bracket_individual_sv(self, vcf_record):
        # Transform REF/ALT to equivalent notation so that REF contains the lowest position
        if vcf_record.alt_sv_bracket.contig == vcf_record.contig and vcf_record.alt_sv_bracket.pos < vcf_record.pos:
            new_contig = vcf_record.alt_sv_bracket.contig
            new_pos = vcf_record.alt_sv_bracket.pos
            new_end = new_pos
            alt_prefix = vcf_record.alt_sv_bracket.suffix
            alt_suffix = vcf_record.alt_sv_bracket.prefix
            alt_bracket = ']' if vcf_record.alt_sv_bracket.bracket == '[' else '['
            alt_contig = vcf_record.contig
            alt_pos = vcf_record.pos
            new_alts = [f'{alt_prefix}{alt_bracket}{alt_contig}:{alt_pos}{alt_bracket}{alt_suffix}']
            alt_sv_bracket = BracketSVRecord(alt_prefix, alt_bracket, alt_contig, alt_pos, alt_suffix)
            vcf_record = VariantRecord(new_contig, new_pos, new_end, vcf_record.id, vcf_record.ref, new_alts,
                                       vcf_record.filter, vcf_record.info, alt_sv_bracket, None)

        # INV -> 1 10 N]1:20] or 1 20 N]1:10]
        #        1 10 [1:20[N or 1 20 [1:10[N
        # DEL -> 1 10 N[1:20[ or 1 20 ]1:10]N
        # DUP -> 1 10 ]1:20]N or 1 20 N[1:10[
        # INS or BND -> any posibility with different contigs
        if vcf_record.contig != vcf_record.alt_sv_bracket.contig:
            # BND & INS with different contig ~ TRN
            return self.__variants.append((VariantType.TRN, vcf_record))
        elif vcf_record.alt_sv_bracket.prefix and vcf_record.alt_sv_bracket.bracket == '[':
            # DEL
            return self.__variants.append((VariantType.DEL, vcf_record))
        elif not vcf_record.alt_sv_bracket.prefix and vcf_record.alt_sv_bracket.bracket == ']':
            # DUP
            return self.__variants.append((VariantType.DUP, vcf_record))
        else:
            # INV
            return self.__variants.append((VariantType.INV, vcf_record))

    def __parse_shorthand_sv(self, vcf_record):
        if vcf_record.alt_sv_shorthand.type == 'DEL':
            return self.__variants.append((VariantType.DEL, vcf_record))
        elif vcf_record.alt_sv_shorthand.type == 'DUP':
            return self.__variants.append((VariantType.DUP, vcf_record))
        elif vcf_record.alt_sv_shorthand.type == 'INV':
            return self.__variants.append((VariantType.INV, vcf_record))
        elif vcf_record.alt_sv_shorthand.type == 'INS':
            return self.__variants.append((VariantType.INS, vcf_record))
        elif vcf_record.alt_sv_shorthand.type == 'CNV':
            return self.__variants.append((VariantType.CNV, vcf_record))

    def __handle_uncertain_sv(self):
        pairs = {}
        paired_records = []
        for alt_name, alt_dicts in self.__pending_sv_pairs.items():
            for sv_name, vcf_record in alt_dicts.items():
                # Check if is already in the dictionary
                sv_id = vcf_record.id
                if sv_id in pairs:
                    # Retrieve saved data
                    previous_record, previous_alt, previous_sv = pairs[sv_id]
                    # Mark for deletion
                    paired_records.append((alt_name, sv_name))
                    paired_records.append((previous_alt, previous_sv))
                    record = _select_record(vcf_record, previous_record)
                    self.__parse_bracket_individual_sv(record)
                    continue

                # Check if labeled with MATEID or PARID
                if 'MATEID' in vcf_record.info:
                    mate_id = vcf_record.info['MATEID']
                elif 'PARID' in vcf_record.info:
                    mate_id = vcf_record.info['PARID']
                mate_id = mate_id[0] if type(mate_id) != str else mate_id
                pairs[mate_id] = (vcf_record, alt_name, sv_name)

        # Remove paired records from pending
        for alt_name, sv_name in paired_records:
            self.__pending_sv_pairs[alt_name].pop(sv_name)
            if len(self.__pending_sv_pairs[alt_name]) == 0:
                self.__pending_sv_pairs.pop(alt_name)
