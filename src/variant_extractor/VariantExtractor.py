# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
from os import path
from argparse import ArgumentParser
import warnings
import pysam

from .private._utils import compare_contigs, permute_bracket_sv, convert_inv_to_bracket
from .private._parser import parse_bracket_sv, parse_shorthand_sv, parse_sgl_sv, parse_standard_record
from .private._PendingBreakends import PendingBreakends
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
        self.__pending_breakends = PendingBreakends()
        # Read the file
        with open(file=vcf_file, mode='r') as vcf_handle:
            save = pysam.set_verbosity(0)
            with pysam.VariantFile(vcf_handle) as vcf:
                self.__variant_file = vcf
                pysam.set_verbosity(save)
                for rec in vcf:
                    self.__handle_record(rec)
        # Handle imprecise unpaired breakends
        self.__handle_imprecise_sv()
        # Only single-paired records or not ensuring pairs
        if not self.ensure_pairs or self.__pairs_found == 0:
            for vcf_record in self.__pending_breakends.values():
                self.__handle_bracket_individual_sv(vcf_record)
        # Found unpaired records
        elif len(self.__pending_breakends) > 0:
            exception_text = ''
            for vcf_record in self.__pending_breakends.values():
                exception_text += f'REF={vcf_record.contig}:{vcf_record.pos}\tALT={vcf_record.alt_sv_bracket.contig}:{vcf_record.alt_sv_bracket.pos}\n'
            raise Exception(
                (f'Unpaired SV breakends:\n{exception_text}'
                 f'Exception: There are {len(self.__pending_breakends)} unpaired SV breakends. '
                 'Please, check the entires show above in VCF file. '
                 'Use ensure_pairs=False to ignore unpaired SV breakends.'))
        return self.__variants

    def __handle_record(self, rec):
        if self.only_pass and 'PASS' not in rec.filter:
            return
        # Handle multiallelic records
        if len(rec.alts) != 1:
            return self.__handle_multiallelic_record(rec)
        # Check if bracket SV record
        vcf_record = parse_bracket_sv(rec)
        if vcf_record:
            return self.__handle_bracket_sv(vcf_record)
        # Check if shorthand SV record
        vcf_record = parse_shorthand_sv(rec)
        if vcf_record:
            return self.__handle_shorthand_sv(vcf_record)
        # Check if single breakend SV record
        vcf_record = parse_sgl_sv(rec)
        if vcf_record:
            return self.__variants.append(vcf_record)
        # Check if standard record
        vcf_record = parse_standard_record(rec)
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
        # Check for pending breakends
        previous_record = self.__pending_breakends.pop(vcf_record)
        if previous_record is None:
            self.__pending_breakends.push(vcf_record)
            return
        # Mate breakend found, handle it
        self.__pairs_found += 1
        if compare_contigs(previous_record.contig, vcf_record.contig) == -1:
            self.__handle_bracket_individual_sv(previous_record)
        else:
            self.__handle_bracket_individual_sv(vcf_record)

    def __handle_bracket_individual_sv(self, vcf_record):
        contig_comparison = compare_contigs(vcf_record.contig, vcf_record.alt_sv_bracket.contig)
        # Transform REF/ALT to equivalent notation so that REF contains the lowest contig and position
        if contig_comparison == 1 or (contig_comparison == 0 and vcf_record.pos > vcf_record.alt_sv_bracket.pos):
            vcf_record = permute_bracket_sv(vcf_record)
        return self.__variants.append(vcf_record)

    def __handle_shorthand_sv(self, vcf_record):
        if vcf_record.variant_type == VariantType.INV:
            # Transform INV into bracket notation
            vcf_record_1, vcf_record_2 = convert_inv_to_bracket(vcf_record)
            self.__variants.append(vcf_record_1)
            self.__variants.append(vcf_record_2)
        else:
            self.__variants.append(vcf_record)

    def __handle_imprecise_sv(self):
        pairs = {}
        paired_records = []
        for vcf_record in self.__pending_breakends.values():
            # Check if is already in the dictionary
            previous_record = pairs.get(vcf_record.id)
            if previous_record:
                del pairs[vcf_record.id]
                # Mark as paired for deletion
                paired_records.append(vcf_record)
                paired_records.append(previous_record)
                record = select_record(vcf_record, previous_record)
                self.__handle_bracket_individual_sv(record)
                continue

            # Check if it has MATEID or PARID
            if 'MATEID' in vcf_record.info:
                mate_id = vcf_record.info['MATEID']
            elif 'PARID' in vcf_record.info:
                mate_id = vcf_record.info['PARID']
            else:
                continue
            # Store record based on mate id
            mate_id = mate_id[0] if type(mate_id) != str else mate_id
            pairs[mate_id] = vcf_record
        # Remove paired records from pending
        for vcf_record in paired_records:
            self.__pending_breakends.remove(vcf_record)

    def __handle_multiallelic_record(self, rec):
        for alt in rec.alts:
            # WARNING: This overrides the record
            rec.alts = [alt]
            self.__handle_record(rec)
