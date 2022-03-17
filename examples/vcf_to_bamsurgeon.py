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

from VariantExtractor import VariantExtractor, VariationType

if __name__ == '__main__':
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
        # TODO: What with INDELs different REF/ALT sizes
        # TODO: Check if 0-based or 1-based for SVs, but probably 0-based
        zero_based_pos = variant_record.pos - 1

        if var_type == VariationType.SNV:
            output_file_snv.write(f'{variant_record.contig} {variant_record.pos} {variant_record.pos} 1.0 {variant_record.alts[0]}\n')
        elif var_type == VariationType.INDEL_DEL:
            output_file_indel.write(
                f'{variant_record.contig} {zero_based_pos} {variant_record.end-1} 1.0 DEL\n')
        elif var_type == VariationType.INDEL_INS:
            output_file_indel.write(
                f'{variant_record.contig} {zero_based_pos} {zero_based_pos+1} 1.0 INS {variant_record.alts[0]}\n')
        else:
            # SVs
            if variant_record.alt_sv_precise:
                alt_contig = variant_record.alt_sv_precise.contig
                zero_based_end_pos = variant_record.alt_sv_precise.pos - 1
            else:
                alt_contig = variant_record.contig
                zero_based_end_pos = variant_record.end - 1
            if var_type == VariationType.TRN:
                # TODO: ++, -- ??
                output_file_sv.write(
                    f'{variant_record.contig} {zero_based_pos} {zero_based_pos} TRN {alt_contig} {zero_based_end_pos} {zero_based_end_pos} ++ 1.0\n')
            elif var_type == VariationType.INV:
                output_file_sv.write(f'{variant_record.contig} {zero_based_pos} {zero_based_end_pos} INV 1.0\n')
            elif var_type == VariationType.DUP:
                output_file_sv.write(f'{variant_record.contig} {zero_based_pos} {zero_based_end_pos} DUP 1.0\n')
            elif var_type == VariationType.DEL:
                output_file_sv.write(f'{variant_record.contig} {zero_based_pos} {zero_based_end_pos} DEL 1.0\n')

    print(f'Reading VCF file: {args.vcf_file}')
    extractor = VariantExtractor(indel_threshold=100)
    for var_type, variant_record in extractor.read_vcf(args.vcf_file):
        variant_callback(var_type, variant_record)


    output_file_sv.close()
    output_file_snv.close()
    output_file_indel.close()
    print(f'Output files generated: {args.output_file_schema}_sv.in, {args.output_file_schema}_snv.in, {args.output_file_schema}_indel.in')
