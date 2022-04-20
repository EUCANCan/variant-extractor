# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
from os import path
from argparse import ArgumentParser
import warnings
import pysam

from ._common import _select_record, _permute_bracket_sv, _convert_inv_to_bracket
from ._parser import _parse_bracket_sv, _parse_shorthand_sv, _parse_sgl_sv, _parse_standard_record
from .variants import VariantType


class VariantExtractor:
    """
    Reads and extracts variants from VCF files. This class is designed to be
    used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis.
    """

    def __init__(self, only_pass=False, ensure_pairs=True):
        """
        Parameters
        ----------
        only_pass : bool, optional
            If :code:`True`, only records with PASS filter will be considered.
        ensure_pairs : bool, optional
            If :code:`True`, throws an exception if a breakend is missing a pair when all other were paired successfully.
        """
        self.ensure_pairs = ensure_pairs
        self.only_pass = only_pass

    def read_vcf(self, vcf_file):
        """Reads VCF file and extracts all variants.

        Parameters
        ----------
        vcf_file : str
            A VCF formatted file. The file is automatically opened.

        Returns
        -------
        list[VariantRecord]
        """

        self.__variants = []
        self.__pairs_found = 0
        self.__pending_sv_pairs = {}
        # Read the file
        with open(file=vcf_file, mode='r') as vcf_handle:
            save = pysam.set_verbosity(0)
            with pysam.VariantFile(vcf_handle) as vcf:
                self.__variant_file = vcf
                pysam.set_verbosity(save)
                for rec in vcf:
                    self.__parse_record(rec)
        # Handle imprecise unpaired SV
        self.__handle_imprecise_sv()
        # Only single-paired records or not ensuring pairs
        if not self.ensure_pairs or self.__pairs_found == 0:
            for alt_dicts in self.__pending_sv_pairs.values():
                for vcf_record in alt_dicts.values():
                    self.__handle_bracket_individual_sv(vcf_record)
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
        if self.only_pass and 'PASS' not in rec.filter:
            return
        # Handle multiallelic records
        if len(rec.alts) != 1:
            return self.__handle_multiallelic_record(rec)
        vcf_record = _parse_bracket_sv(rec)
        # Check if bracket SV record
        if vcf_record:
            return self.__handle_bracket_sv(vcf_record)
        # Check if shorthand SV record
        vcf_record = _parse_shorthand_sv(rec)
        if vcf_record:
            return self.__handle_shorthand_sv(vcf_record)
        # Check if single breakend SV record
        vcf_record = _parse_sgl_sv(rec)
        if vcf_record:
            return self.__variants.append(vcf_record)
        # Check if standard record
        vcf_record = _parse_standard_record(rec)
        if vcf_record:
            return self.__handle_standard_record(vcf_record)
        else:
            warnings.warn(f'Skipping unrecognized record:\n{rec}')

    def __handle_standard_record(self, vcf_record):
        # Normalize complex indels
        i = 0
        min_size = min(len(vcf_record.ref), len(vcf_record.alt)) - 1
        while i < min_size:
            # Atomize SNVs
            if vcf_record.ref[i] != vcf_record.alt[i]:
                new_vcf_record = vcf_record._replace(
                    ref=vcf_record.ref[i], pos=i+vcf_record.pos, end=i+vcf_record.pos, length=0, alt=vcf_record.alt[i], id=f'{vcf_record.id}_{i}' if vcf_record.id else None, variant_type=VariantType.SNV)
                self.__variants.append(new_vcf_record)
            i += 1

        if len(vcf_record.ref) > len(vcf_record.alt):
            # Deletion
            new_vcf_record = vcf_record._replace(
                ref=vcf_record.ref[i:], pos=i+vcf_record.pos, end=i+vcf_record.pos, length=len(vcf_record.ref)-i-1, alt=vcf_record.ref[i], id=f'{vcf_record.id}_{i}' if vcf_record.id else None, variant_type=VariantType.DEL)
            self.__variants.append(new_vcf_record)
        elif len(vcf_record.ref) < len(vcf_record.alt):
            # Insertion
            new_vcf_record = vcf_record._replace(
                ref=vcf_record.ref[i:], pos=i+vcf_record.pos, end=i+vcf_record.pos, length=len(vcf_record.alt)-i-1, alt=vcf_record.ref[i]+vcf_record.alt[i+1:], id=f'{vcf_record.id}_{i}' if vcf_record.id else None, variant_type=VariantType.INS)
            self.__variants.append(new_vcf_record)

        if vcf_record.ref[i] != vcf_record.alt[i]:
            # Setup record id
            id_offset = 0 if len(vcf_record.ref) == len(vcf_record.alt) else 1
            id_str = f'_{i+id_offset}' if i > 0 else ''
            new_record_id = vcf_record.id+id_str if vcf_record.id else None
            # Last SNV
            new_vcf_record = vcf_record._replace(
                ref=vcf_record.ref[i], pos=i+vcf_record.pos, end=i+vcf_record.pos, length=0, alt=vcf_record.alt[i], id=new_record_id, variant_type=VariantType.SNV)
            self.__variants.append(new_vcf_record)

    def __handle_bracket_sv(self, vcf_record):
        # Check for pending SVs
        previous_record = self.__pop_pending_sv_pair(vcf_record)
        if previous_record is None:
            self.__store_pending_sv_pair(vcf_record)
            return
        # Mate SV found, parse it
        self.__pairs_found += 1
        record = _select_record(previous_record, vcf_record)
        self.__handle_bracket_individual_sv(record)

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

    def __handle_bracket_individual_sv(self, vcf_record):
        # Transform REF/ALT to equivalent notation so that REF contains the lowest position
        if vcf_record.alt_sv_bracket.contig == vcf_record.contig and vcf_record.alt_sv_bracket.pos < vcf_record.pos:
            vcf_record = _permute_bracket_sv(vcf_record)
        return self.__variants.append(vcf_record)

    def __handle_shorthand_sv(self, vcf_record):
        if vcf_record.variant_type == VariantType.INV:
            # Transform INV into bracket notation
            vcf_record_1, vcf_record_2 = _convert_inv_to_bracket(vcf_record)
            self.__variants.append(vcf_record_1)
            self.__variants.append(vcf_record_2)
        else:
            self.__variants.append(vcf_record)

    def __handle_imprecise_sv(self):
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
                    self.__handle_bracket_individual_sv(record)
                    continue

                # Check if labeled with MATEID or PARID
                if 'MATEID' in vcf_record.info:
                    mate_id = vcf_record.info['MATEID']
                elif 'PARID' in vcf_record.info:
                    mate_id = vcf_record.info['PARID']
                else:
                    continue
                mate_id = mate_id[0] if type(mate_id) != str else mate_id
                pairs[mate_id] = (vcf_record, alt_name, sv_name)

        # Remove paired records from pending
        for alt_name, sv_name in paired_records:
            self.__pending_sv_pairs[alt_name].pop(sv_name)
            if len(self.__pending_sv_pairs[alt_name]) == 0:
                self.__pending_sv_pairs.pop(alt_name)

    def __handle_multiallelic_record(self, rec):
        for alt in rec.alts:
            # WARNING: This overrides the record
            rec.alts = [alt]
            self.__parse_record(rec)
