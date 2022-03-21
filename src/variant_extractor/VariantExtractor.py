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

# TODO: What to do with REF=N ALT=.N
# Regex for SVs
PRECISE_SV_REGEX = re.compile(r'([A-Za-z]*)(\[|\])([^\]\[:]+:[0-9]+)(\[|\])([A-Za-z]*)')
NAMED_SV_REGEX = re.compile(r'<(DEL|INS|DUP|INV|CNV])(:[A-Z]+)*>')
NUMBER_CONTIG_REGEX = re.compile(r'[0-9]+')


class VariationType(Enum):
    SNV = auto()
    INDEL_INS = auto()
    INDEL_DEL = auto()
    DEL = auto()
    INS = auto()
    DUP = auto()
    INV = auto()
    CNV = auto()
    TRN = auto()


class PreciseSVRecord(NamedTuple):
    prefix: Optional[str]
    bracket: str
    contig: str
    pos: int
    suffix: Optional[str]


class NamedSVRecord(NamedTuple):
    type: str
    extra: List[str]


class VariantRecord(NamedTuple):
    contig: str
    pos: int
    end: int
    id: str
    ref: str
    alts: List[str]
    filter: str
    info: dict
    alt_sv_precise: Optional[PreciseSVRecord]  # G]17:198982]
    alt_sv_named: Optional[NamedSVRecord]  # <DEL:ME:ALU>


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


class VariantExtractor:
    def __init__(self, indel_threshold=-1, ensure_pairs=True):
        self.indel_threshold = indel_threshold
        self.ensure_pairs = ensure_pairs

    def read_vcf(self, vcf_file):
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
                    self.__parse_precise_individual_sv(vcf_record)
        # Found unpaired records
        elif len(self.__pending_sv_pairs) > 0:
            exception_text = ''
            for alt_dicts in self.__pending_sv_pairs.values():
                for vcf_record in alt_dicts.values():
                    exception_text += f'REF={vcf_record.contig}:{vcf_record.pos}\tALT={vcf_record.alt_sv_precise.contig}:{vcf_record.alt_sv_precise.pos}\n'
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
        vcf_record = VariantRecord(rec.contig, rec.pos, rec.stop, rec.id, rec.ref,
                                   rec.alts, rec.filter, rec.info, None, None)
        sv_match_precise = PRECISE_SV_REGEX.search(rec.alts[0])
        sv_match_named = NAMED_SV_REGEX.search(rec.alts[0])
        # Check if precise SV
        if sv_match_precise:
            # Extract data from regex
            alt_prefix = sv_match_precise.group(1)
            alt_bracket = sv_match_precise.group(2)
            alt_contig, alt_pos = sv_match_precise.group(3).split(':')
            alt_suffix = sv_match_precise.group(5)
            alt_sv_precise = PreciseSVRecord(alt_prefix, alt_bracket, alt_contig, int(alt_pos), alt_suffix)
            vcf_record = vcf_record._replace(alt_sv_precise=alt_sv_precise)
            return self.__parse_precise_sv(vcf_record)
        # Check if named SV
        elif sv_match_named:
            # Extract data from regex
            alt_type = sv_match_named.group(1)
            alt_extra = sv_match_named.group(2).split(':') if sv_match_named.group(2) else None
            alt_sv_named = NamedSVRecord(alt_type, alt_extra)
            vcf_record = vcf_record._replace(alt_sv_named=alt_sv_named)
            return self.__parse_named_sv(vcf_record)
        # Check if SNV
        if len(vcf_record.alts[0]) == len(vcf_record.ref):
            # REF=CTT ALT=ATG -> Normalize to 3 SNVs
            for i in range(len(vcf_record.ref)):
                if vcf_record.alts[0][i] != vcf_record.ref[i]:
                    new_vcf_record = vcf_record._replace(
                        ref=vcf_record.ref[i], pos=i+vcf_record.pos, end=i+vcf_record.end, alts=[vcf_record.alts[0][i]])
                    self.__variants.append((VariationType.SNV, new_vcf_record))
        # Check if INDEL_DEL
        elif len(vcf_record.alts[0]) < len(vcf_record.ref):
            if self.indel_threshold != -1 and abs(vcf_record.pos - vcf_record.end) > self.indel_threshold:
                return self.__variants.append((VariationType.DEL, vcf_record))
            else:
                return self.__variants.append((VariationType.INDEL_DEL, vcf_record))
        # Check if INDEL_INS
        else:
            if self.indel_threshold != -1 and len(vcf_record.alts[0]) > self.indel_threshold:
                return self.__variants.append((VariationType.INS, vcf_record))
            else:
                return self.__variants.append((VariationType.INDEL_INS, vcf_record))

    def __parse_precise_sv(self, vcf_record):
        # Check for pending SVs
        previous_record = self.__pop_pending_sv_pair(vcf_record)
        if previous_record is None:
            self.__store_pending_sv_pair(vcf_record)
            return
        # Mate SV found, parse it
        self.__pairs_found += 1
        record = select_record(previous_record, vcf_record)
        self.__parse_precise_individual_sv(record)

    def __store_pending_sv_pair(self, vcf_record):
        alt_name = f'{vcf_record.alt_sv_precise.contig}{vcf_record.alt_sv_precise.pos}'
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
        previous_record_alt_name = f'{vcf_record.alt_sv_precise.contig}{vcf_record.alt_sv_precise.pos}'
        previous_record_alt = previous_records.get(previous_record_alt_name)
        if previous_record_alt is None:
            return None
        previous_records.pop(previous_record_alt_name)
        if len(previous_records) == 0:
            self.__pending_sv_pairs.pop(sv_name)
        return previous_record_alt

    def __parse_precise_individual_sv(self, vcf_record):
        # Transform REF/ALT to equivalent notation so that REF contains the lowest position
        if vcf_record.alt_sv_precise.contig == vcf_record.contig and vcf_record.alt_sv_precise.pos < vcf_record.pos:
            new_contig = vcf_record.alt_sv_precise.contig
            new_pos = vcf_record.alt_sv_precise.pos
            new_end = new_pos
            alt_prefix = vcf_record.alt_sv_precise.suffix
            alt_suffix = vcf_record.alt_sv_precise.prefix
            alt_bracket = ']' if vcf_record.alt_sv_precise.bracket == '[' else '['
            alt_contig = vcf_record.contig
            alt_pos = vcf_record.pos
            new_alts = [f'{alt_prefix}{alt_bracket}{alt_contig}:{alt_pos}{alt_bracket}{alt_suffix}']
            alt_sv_precise = PreciseSVRecord(alt_prefix, alt_bracket, alt_contig, alt_pos, alt_suffix)
            vcf_record = VariantRecord(new_contig, new_pos, new_end, vcf_record.id, vcf_record.ref, new_alts,
                                   vcf_record.filter, vcf_record.info, alt_sv_precise, None)

        # TODO: Edit end property for DEL, INV and DUP

        # INV -> 1 10 N]1:20] or 1 20 N]1:10]
        #        1 10 [1:20[N or 1 20 [1:10[N
        # DEL -> 1 10 N[1:20[ or 1 20 ]1:10]N
        # DUP -> 1 10 ]1:20]N or 1 20 N[1:10[
        # INS or BND -> any posibility with different contigs
        if vcf_record.contig != vcf_record.alt_sv_precise.contig:
            # BND & INS with different contig ~ TRN
            return self.__variants.append((VariationType.TRN, vcf_record))
        elif vcf_record.alt_sv_precise.prefix and vcf_record.alt_sv_precise.bracket == '[':
            # DEL
            if self.indel_threshold != -1 and abs(vcf_record.alt_sv_precise.pos - vcf_record.pos) < self.indel_threshold:
                return self.__variants.append((VariationType.INDEL_DEL, vcf_record))
            else:
                return self.__variants.append((VariationType.DEL, vcf_record))
        elif not vcf_record.alt_sv_precise.prefix and vcf_record.alt_sv_precise.bracket == ']':
            # DUP
            return self.__variants.append((VariationType.DUP, vcf_record))
        else:
            # INV
            return self.__variants.append((VariationType.INV, vcf_record))

    def __parse_named_sv(self, vcf_record):
        if vcf_record.alt_sv_named.type == 'DEL':
            return self.__variants.append((VariationType.DEL, vcf_record))
        elif vcf_record.alt_sv_named.type == 'DUP':
            return self.__variants.append((VariationType.DUP, vcf_record))
        elif vcf_record.alt_sv_named.type == 'INV':
            return self.__variants.append((VariationType.INV, vcf_record))
        elif vcf_record.alt_sv_named.type == 'INS':
            return self.__variants.append((VariationType.INS, vcf_record))
        elif vcf_record.alt_sv_named.type == 'CNV':
            return self.__variants.append((VariationType.CNV, vcf_record))

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
                    record = select_record(vcf_record, previous_record)
                    self.__parse_precise_individual_sv(record)
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
