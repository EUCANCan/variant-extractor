# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# MIT License
'''
Generates the BAMSurgeon input from a VCF file
Expected usage:
    $ python vcf_to_bamsurgeon.py <vcf_file> <output_file_schema>
Use --help for more information.
'''
from os import path
from argparse import ArgumentParser
import re

VAF = 0.5
INDEL_THRESHOLD = 100

if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + '/../src/')
    from variant_extractor import VariantExtractor, VariantType

    # Parse arguments
    parser = ArgumentParser(description='Generate BAMSurgeon input from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_file_schema', help='Output file schema')
    args = parser.parse_args()

    # Open output files
    output_file_sv = open(f'{args.output_file_schema}_sv.in', 'w')
    output_file_snv = open(f'{args.output_file_schema}_snv.in', 'w')
    output_file_indel = open(f'{args.output_file_schema}_indel.in', 'w')

    def variant_callback(var_type, variant_record):
        if var_type == VariantType.SNV:
            output_file_snv.write(
                f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {VAF} {variant_record.alts[0]}\n')
        elif var_type == VariantType.INDEL_DEL:
            # TODO: What with indels different REF/ALT sizes
            # Must be 0-based
            zero_based_pos = variant_record.pos - 1
            zero_based_end = variant_record.end - 1
            output_file_indel.write(
                f'{variant_record.contig} {zero_based_pos} {zero_based_end} {VAF} DEL\n')
        elif var_type == VariantType.INDEL_INS:
            pass
            # TODO: What to insert if not defined
            # zero_based_pos = variant_record.pos - 1
            # output_file_indel.write(
            #     f'{variant_record.contig} {zero_based_pos} {zero_based_pos+1} {VAF} INS {variant_record.alts[0]}\n')
        else:
            if var_type == VariantType.TRN:
                alt_contig = variant_record.alt_sv_bracket.contig
                alt_pos = variant_record.alt_sv_bracket.pos
                # Calculate strand notation
                if variant_record.alt_sv_bracket.bracket == '[' and variant_record.alt_sv_bracket.prefix:
                    strand_notation = '++'
                elif variant_record.alt_sv_bracket.bracket == ']' and variant_record.alt_sv_bracket.prefix:
                    strand_notation = '+-'
                elif variant_record.alt_sv_bracket.bracket == '[' and not variant_record.alt_sv_bracket.prefix:
                    strand_notation = '-+'
                else:
                    strand_notation = '--'

                output_file_sv.write(
                    f'{variant_record.contig} {variant_record.pos} {variant_record.pos} TRN {alt_contig} {alt_pos} {alt_pos} {strand_notation} {VAF}\n')
            elif var_type == VariantType.INV:
                output_file_sv.write(f'{variant_record.contig} {variant_record.pos} {variant_record.end} INV {VAF}\n')
            elif var_type == VariantType.DUP:
                output_file_sv.write(f'{variant_record.contig} {variant_record.pos} {variant_record.end} DUP 1 {VAF}\n')
            elif var_type == VariantType.DEL:
                output_file_sv.write(f'{variant_record.contig} {variant_record.pos} {variant_record.end} DEL {VAF}\n')
                # TODO: Add indel threshold
            elif var_type == VariantType.INS:
                # TODO: What to insert if not defined
                pass
                # output_file_sv.write(f'{variant_record.contig} {variant_record.pos} {variant_record.pos} INS {VAF}\n')

    print(f'Reading VCF file: {args.vcf_file}')
    extractor = VariantExtractor()
    for var_type, variant_record in extractor.read_vcf(args.vcf_file):
        variant_callback(var_type, variant_record)

    output_file_sv.close()
    output_file_snv.close()
    output_file_indel.close()
    print(
        f'Output files generated: {args.output_file_schema}_sv.in, {args.output_file_schema}_snv.in, {args.output_file_schema}_indel.in')
