# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# MIT License
from typing import List, Optional
import warnings
import pysam

from .private._utils import compare_contigs, permute_breakend_sv, convert_inv_to_breakend, convert_del_to_ins
from .private._parser import parse_breakend_sv, parse_shorthand_sv, parse_sgl_sv, parse_standard_record
from .private._PendingBreakends import PendingBreakends
from .variants import VariantType
from .variants import VariantRecord

DATAFRAME_COLUMNS = ['start_chrom', 'start', 'end_chrom', 'end', 'ref',
                     'alt', 'length', 'brackets', 'type_inferred', 'variant_record_obj']


class VariantExtractor:
    """
    Reads and extracts variants from VCF files. This class is designed to be
    used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis.
    """

    @staticmethod
    def empty_dataframe():
        """Returns an empty pandas DataFrame with the columns used by this class.
        """
        import pandas as pd
        return pd.DataFrame(columns=DATAFRAME_COLUMNS)

    def __init__(self, vcf_file: str, pass_only=False, ensure_pairs=True, fasta_ref: Optional[str] = None):
        """
        Parameters
        ----------
        vcf_file : str
            A VCF formatted file. The file is automatically opened.
        pass_only : bool, optional
            If :code:`True`, only records with PASS filter will be considered.
        ensure_pairs : bool, optional
            If :code:`True`, throws an exception if a breakend is missing a pair when all other were paired successfully.
        fasta_ref : str, optional
            A FASTA file with the reference genome. Must be indexed.
        """
        self.__ensure_pairs = ensure_pairs
        self.__pass_only = pass_only
        self.__pairs_found = 0
        self.__pending_breakends = PendingBreakends()
        self.__fasta_ref = None
        # Open FASTA file
        if fasta_ref is not None:
            self.__fasta_ref = pysam.FastaFile(fasta_ref)
        # Open VCF file
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
        # Remove non-PASS records from the pending breakends if pass_only is True
        if self.__pass_only:
            vcf_records = list(self.__pending_breakends.values())
            for vcf_record in vcf_records:
                if 'PASS' not in vcf_record.filter:
                    self.__pending_breakends.remove(vcf_record)
        # Only single-paired records or not ensuring pairs
        if not self.__ensure_pairs or self.__pairs_found == 0:
            for vcf_record in self.__pending_breakends.values():
                yield from self.__handle_breakend_individual_sv(vcf_record)
        # Found unpaired records
        elif len(self.__pending_breakends) > 0:
            exception_text = ''
            for vcf_record in self.__pending_breakends.values():
                exception_text += str(vcf_record)+'\n'
            raise Exception(
                (f'There are {len(self.__pending_breakends)} unpaired SV breakends. '
                 'Please, check the entires shown below in VCF file. '
                 f'Use ensure_pairs=False to ignore unpaired SV breakends.\n{exception_text}'))

    def __handle_record(self, rec: pysam.VariantRecord) -> List[VariantRecord]:
        if not rec.alts:
            return []
        if not rec.ref:
            raise ValueError('Record does not have a REF field')
        # Handle multiallelic records
        if len(rec.alts) != 1:
            return self.__handle_multiallelic_record(rec)
        # Check if breakend SV record
        vcf_record = parse_breakend_sv(rec)
        if vcf_record:
            return self.__handle_breakend_sv(vcf_record)
        # Check PASS filter
        if self.__pass_only and 'PASS' not in rec.filter:
            return []
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
                        ref=vcf_record.ref[i], pos=i+vcf_record.pos, end=i+vcf_record.pos, length=1, alt=vcf_record.alt[i], id=variant_id, variant_type=VariantType.SNV)
                    record_list.append(new_vcf_record)
        elif len(vcf_record.ref) > len(vcf_record.alt):
            # Deletion
            vcf_record.variant_type=VariantType.DEL
            record_list.append(vcf_record)
        elif len(vcf_record.ref) < len(vcf_record.alt):
            # Insertion
            vcf_record.variant_type=VariantType.INS
            record_list.append(vcf_record)
        return record_list

    def __handle_breakend_sv(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        # Check for pending breakends
        previous_record = self.__pending_breakends.pop(vcf_record)
        if previous_record is None:
            self.__pending_breakends.push(vcf_record)
            return []
        # Mate breakend found, handle it
        self.__pairs_found += 1
        return self.__handle_braked_paired_sv(previous_record, vcf_record)

    def __handle_braked_paired_sv(self, vcf_record_1: VariantRecord, vcf_record_2: VariantRecord) -> List[VariantRecord]:
        # Check PASS filter
        if self.__pass_only and ('PASS' not in vcf_record_1.filter or 'PASS' not in vcf_record_2.filter):
            return []
        # Unify filters
        filters = set(vcf_record_1.filter) | set(vcf_record_2.filter)
        filters.discard('PASS')
        if len(filters) > 0:
            vcf_record_1.filter = list(filters)
            vcf_record_2.filter = list(filters)
        contig_comparison = compare_contigs(vcf_record_1.contig, vcf_record_2.contig)
        if contig_comparison == 0:
            if vcf_record_1.pos < vcf_record_2.pos:
                return self.__handle_breakend_individual_sv(vcf_record_1)
            else:
                return self.__handle_breakend_individual_sv(vcf_record_2)
        elif contig_comparison == -1:
            return self.__handle_breakend_individual_sv(vcf_record_1)
        else:
            return self.__handle_breakend_individual_sv(vcf_record_2)

    def __handle_breakend_individual_sv(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        assert vcf_record.alt_sv_breakend is not None
        contig_comparison = compare_contigs(vcf_record.contig, vcf_record.alt_sv_breakend.contig)
        # Transform REF/ALT to equivalent notation so that REF contains the lowest contig and position
        if contig_comparison == 1 or (contig_comparison == 0 and vcf_record.pos > vcf_record.alt_sv_breakend.pos):
            vcf_record = permute_breakend_sv(vcf_record, self.__fasta_ref)
        # Handle DEL notated variant as INS
        if vcf_record.length == 1 and vcf_record.variant_type == VariantType.DEL and vcf_record.alt_sv_breakend is not None:
            vcf_record = convert_del_to_ins(vcf_record)
            if vcf_record.length == 0:
                return []
        return [vcf_record]

    def __handle_shorthand_sv(self, vcf_record: VariantRecord) -> List[VariantRecord]:
        if vcf_record.variant_type == VariantType.INV:
            # Transform INV into breakend notation
            vcf_record_1, vcf_record_2 = convert_inv_to_breakend(vcf_record, self.__fasta_ref)
            return [vcf_record_1, vcf_record_2]
        else:
            return [vcf_record]

    def __handle_multiallelic_record(self, rec: pysam.VariantRecord) -> List[VariantRecord]:
        record_list = []
        fake_rec = rec.copy()
        assert fake_rec.alts is not None and len(fake_rec.alts) > 1
        alts = fake_rec.alts
        samples = dict()
        for sample_name in rec.samples:
            sample_dict = dict()
            for key, value in rec.samples[sample_name].items():
                sample_dict[key] = value
            samples[sample_name] = sample_dict

        original_id = fake_rec.id
        for i, alt in enumerate(alts):
            # WARNING: This overrides the record
            fake_rec.alts = (alt,)
            if original_id:
                new_id = f'{original_id}_{i}'
                fake_rec.id = new_id
            new_records = self.__handle_record(fake_rec)
            new_samples = dict()
            for sample_name in samples:
                new_samples[sample_name] = dict()
                for key, value in samples[sample_name].items():
                    if not hasattr(value, '__iter__') or len(value) == self.__variant_file.header.formats[key].number:
                        new_samples[sample_name][key] = value
                    else:
                        if key == 'GT':
                            new_samples[sample_name][key] = (0, samples[sample_name][key][1])
                        elif len(value) == len(alts) + 1:
                            new_samples[sample_name][key] = (value[0], value[i + 1])
                        elif hasattr(value, '__iter__') and len(value) % len(alts) == 0:
                            new_samples[sample_name][key] = value[i::len(alts)]
                        else:
                            new_samples[sample_name][key] = value
            for new_record in new_records:
                new_record.samples = new_samples
            record_list.extend(new_records)
        return record_list

    def to_dataframe(self):
        import pandas as pd
        variants = []

        for variant_record in self:
            start_chrom = variant_record.contig.replace('chr', '')
            start = variant_record.pos
            ref = variant_record.ref
            alt = variant_record.alt
            length = variant_record.length
            end = variant_record.end
            if variant_record.alt_sv_breakend:
                end_chrom = variant_record.alt_sv_breakend.contig.replace('chr', '')
                if start_chrom != end_chrom:
                    end = variant_record.alt_sv_breakend.pos
            else:
                end_chrom = start_chrom

            # Inferred type
            type_inferred = variant_record.variant_type.name
            # breakends
            breakends = ''
            if type_inferred == VariantType.DEL.name:
                breakends = 'N['
            elif type_inferred == VariantType.DUP.name:
                breakends = ']N'
            elif type_inferred == VariantType.INV.name:
                assert variant_record.alt_sv_breakend is not None
                prefix = 'N' if variant_record.alt_sv_breakend.prefix else ''
                suffix = 'N' if variant_record.alt_sv_breakend.suffix else ''
                breakends = prefix + variant_record.alt_sv_breakend.bracket + suffix
            elif type_inferred == VariantType.TRA.name:
                assert variant_record.alt_sv_breakend is not None
                prefix = 'N' if variant_record.alt_sv_breakend.prefix else ''
                suffix = 'N' if variant_record.alt_sv_breakend.suffix else ''
                breakends = prefix + variant_record.alt_sv_breakend.bracket + variant_record.alt_sv_breakend.bracket + suffix

            variants.append([start_chrom, start, end_chrom, end, ref, alt,
                            length, breakends, type_inferred, variant_record])

        return pd.DataFrame(variants, columns=DATAFRAME_COLUMNS)
