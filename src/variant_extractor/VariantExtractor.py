# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
from typing import Generator, List
from os import path
from argparse import ArgumentParser
import warnings
import pysam

from .private._utils import compare_contigs, permute_bracket_sv, convert_inv_to_bracket
from .private._parser import parse_bracket_sv, parse_shorthand_sv, parse_sgl_sv, parse_standard_record
from .private._PendingBreakends import PendingBreakends
from .variants import VariantType
from .variants import VariantRecord


class VariantExtractor:
    """
    Reads and extracts variants from VCF files. This class is designed to be
    used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis.
    """

    def __init__(self, vcf_file: str, pass_only=False, ensure_pairs=True):
        """
        Parameters
        ----------
        vcf_file : str
            A VCF formatted file. The file is automatically opened.
        pass_only : bool, optional
            If :code:`True`, only records with PASS filter will be considered.
        ensure_pairs : bool, optional
            If :code:`True`, throws an exception if a breakend is missing a pair when all other were paired successfully.
        """
        self.__ensure_pairs = ensure_pairs
        self.__pass_only = pass_only
        self.__pairs_found = 0
        self.__pending_breakends = PendingBreakends()
        # Open the file
        vcf_handle = open(file=vcf_file, mode='r')
        save = pysam.set_verbosity(0)
        self.__variant_file = pysam.VariantFile(vcf_handle)
        pysam.set_verbosity(save)

    def close(self):
        """Closes the VCF file.
        """
        self.__variant_file.close()

    def __iter__(self):
        # Read the next record from the VCF file
        for rec in self.__variant_file:
            yield from self.__handle_record(rec)
        # Handle imprecise unpaired breakends
        yield from self.__handle_imprecise_sv()
        # Only single-paired records or not ensuring pairs
        if not self.__ensure_pairs or self.__pairs_found == 0:
            for vcf_record in self.__pending_breakends.values():
                yield from self.__handle_bracket_individual_sv(vcf_record)
        # Found unpaired records
        elif len(self.__pending_breakends) > 0:
            exception_text = ''
            for vcf_record in self.__pending_breakends.values():
                exception_text += str(vcf_record)+'\n'
            raise Exception(
                (f'Unpaired SV breakends:\n{exception_text}'
                 f'Exception: There are {len(self.__pending_breakends)} unpaired SV breakends. '
                 'Please, check the entires show above in VCF file. '
                 'Use ensure_pairs=False to ignore unpaired SV breakends.'))

    def __handle_record(self, rec: pysam.VariantRecord) -> List[VariantRecord]:
        if self.__pass_only and 'PASS' not in rec.filter:
            return []
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
            return [vcf_record]
        # Check if standard record
        vcf_record = parse_standard_record(rec)
        if vcf_record:
            return self.__handle_standard_record(vcf_record)
        else:
            warnings.warn(f'Skipping unrecognized record:\n{rec}')
            return []

    def __handle_standard_record(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        record_list = []
        if len(vcf_record.ref) == len(vcf_record.alt):
            for i in range(len(vcf_record.ref)):
                # Atomize SNVs
                if vcf_record.ref[i] != vcf_record.alt[i]:
                    variant_id = f'{vcf_record.id}_{i}' if len(vcf_record.ref) > 1 and vcf_record.id is not None \
                        else vcf_record.id
                    new_vcf_record = vcf_record._replace(
                        ref=vcf_record.ref[i], pos=i+vcf_record.pos, end=i+vcf_record.pos, length=0, alt=vcf_record.alt[i], id=variant_id, variant_type=VariantType.SNV)
                    record_list.append(new_vcf_record)
        elif len(vcf_record.ref) > len(vcf_record.alt):
            # Deletion
            new_vcf_record = vcf_record._replace(variant_type=VariantType.DEL)
            record_list.append(new_vcf_record)
        elif len(vcf_record.ref) < len(vcf_record.alt):
            # Insertion
            new_vcf_record = vcf_record._replace(variant_type=VariantType.INS)
            record_list.append(new_vcf_record)
        return record_list

    def __handle_bracket_sv(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        # Check for pending breakends
        previous_record = self.__pending_breakends.pop(vcf_record)
        if previous_record is None:
            self.__pending_breakends.push(vcf_record)
            return []
        # Mate breakend found, handle it
        self.__pairs_found += 1
        contig_comparison = compare_contigs(previous_record.contig, vcf_record.contig)
        if contig_comparison == 0:
            if previous_record.pos < vcf_record.pos:
                return self.__handle_bracket_individual_sv(previous_record)
            else:
                return self.__handle_bracket_individual_sv(vcf_record)
        elif contig_comparison == -1:
            return self.__handle_bracket_individual_sv(previous_record)
        elif contig_comparison == 1:
            return self.__handle_bracket_individual_sv(vcf_record)

    def __handle_bracket_individual_sv(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        contig_comparison = compare_contigs(vcf_record.contig, vcf_record.alt_sv_bracket.contig)
        # Transform REF/ALT to equivalent notation so that REF contains the lowest contig and position
        if contig_comparison == 1 or (contig_comparison == 0 and vcf_record.pos > vcf_record.alt_sv_bracket.pos):
            vcf_record = permute_bracket_sv(vcf_record)
        return [vcf_record]

    def __handle_shorthand_sv(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        if vcf_record.variant_type == VariantType.INV:
            # Transform INV into bracket notation
            vcf_record_1, vcf_record_2 = convert_inv_to_bracket(vcf_record)
            return [vcf_record_1, vcf_record_2]
        else:
            return [vcf_record]

    def __handle_imprecise_sv(self) -> Generator[VariantRecord, None, None]:
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
                if compare_contigs(previous_record.contig, vcf_record.contig) == -1:
                    yield from self.__handle_bracket_individual_sv(previous_record)
                else:
                    yield from self.__handle_bracket_individual_sv(vcf_record)
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

    def __handle_multiallelic_record(self, rec: pysam.VariantRecord) -> List[VariantRecord]:
        record_list = []
        alts = rec.alts
        original_id = rec.id
        for i, alt in enumerate(alts):
            # WARNING: This overrides the record
            rec.alts = [alt]
            new_id = f'{original_id}_{i}'
            rec.id = new_id
            record_list += self.__handle_record(rec)
        return record_list
