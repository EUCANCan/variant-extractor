# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# BSC Dual License
'''
Generates the BAMSurgeon input from a VCF file
Expected usage:
    $ python vcf_to_bamsurgeon.py <vcf_file> <output_file_schema>
Use --help for more information.
'''
from argparse import ArgumentParser
import random

VAF = 0.5
INDEL_THRESHOLD = 90


def _generate_random_dna(length):
    '''
    Generates a random DNA sequence of the given length
    '''
    return ''.join(random.choices(['A', 'C', 'G', 'T'], weights=[0.3, 0.2, 0.2, 0.3], k=length))


if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + '/../src/')
    from variant_extractor import VariantExtractor
    from variant_extractor.variants import VariantType

    # Parse arguments
    parser = ArgumentParser(description='Generate BAMSurgeon input from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_file_schema', help='Output file schema')
    args = parser.parse_args()

    # Open output files
    output_file_sv = open(f'{args.output_file_schema}_sv.in', 'w')
    output_file_snv = open(f'{args.output_file_schema}_snv.in', 'w')
    output_file_indel = open(f'{args.output_file_schema}_indel.in', 'w')
    
    vaf = VAF
    indel_threshold = INDEL_THRESHOLD

    print(f'Reading VCF file: {args.vcf_file}')
    extractor = VariantExtractor(args.vcf_file, pass_only=True)
    for variant_record in extractor:
        if variant_record.variant_type == VariantType.SNV:
            output_file_snv.write(
                f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {vaf} {variant_record.alt}\n')
        else:
            # Add prefix or suffix as insertion. Ex: AAAGGTC[1:12121[
            insertion_prefix = ''
            if variant_record.alt_sv_breakend:
                insertion_prefix = f'INS {variant_record.alt_sv_breakend.prefix[1:]};' if variant_record.alt_sv_breakend.prefix is not None and len(
                    variant_record.alt_sv_breakend.prefix) > 1 else ''
                insertion_prefix = f'INS {variant_record.alt_sv_breakend.suffix[:-1]};' if variant_record.alt_sv_breakend.suffix is not None and len(
                    variant_record.alt_sv_breakend.suffix) > 1 else ''

            if variant_record.variant_type == VariantType.TRN or variant_record.variant_type == VariantType.INV:
                assert variant_record.alt_sv_breakend is not None
                # Convert INV to TRN since most of them are not complete
                alt_contig = variant_record.alt_sv_breakend.contig
                alt_pos = variant_record.alt_sv_breakend.pos
                # Calculate strand notation
                if variant_record.alt_sv_breakend.bracket == '[' and variant_record.alt_sv_breakend.prefix:
                    strand_notation = '++'
                elif variant_record.alt_sv_breakend.bracket == ']' and variant_record.alt_sv_breakend.prefix:
                    strand_notation = '+-'
                elif variant_record.alt_sv_breakend.bracket == '[' and not variant_record.alt_sv_breakend.prefix:
                    strand_notation = '-+'
                else:
                    strand_notation = '--'

                op = f'TRN {alt_contig} {alt_pos} {alt_pos} {strand_notation} {vaf}'
                output_file_sv.write(
                    f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {insertion_prefix}{op}\n')
            elif variant_record.variant_type == VariantType.DUP:
                op = f'DUP 1 {vaf}'
                output_file_sv.write(
                    f'{variant_record.contig} {variant_record.pos} {variant_record.end} {insertion_prefix}{op}\n')
            elif variant_record.variant_type == VariantType.DEL:
                # Check if INDEL
                if variant_record.end - variant_record.pos < indel_threshold:
                    output_file_indel.write(
                        f'{variant_record.contig} {variant_record.pos-1} {variant_record.end-1} {vaf} DEL\n')
                else:
                    op = f'DEL {vaf}'
                    output_file_sv.write(
                        f'{variant_record.contig} {variant_record.pos} {variant_record.end} {insertion_prefix}{op}\n')
            elif variant_record.variant_type == VariantType.INS:
                if variant_record.alt_sv_shorthand:
                    insert_length = variant_record.length
                    dna_sequence = _generate_random_dna(insert_length)
                else:
                    dna_sequence = variant_record.alt
                    insert_length = len(dna_sequence)
                # Check if INDEL
                if insert_length < indel_threshold:
                    output_file_indel.write(
                        f'{variant_record.contig} {variant_record.pos-1} {variant_record.pos} {vaf} INS {dna_sequence}\n')
                else:
                    # Cannot set VAF for SV insertions
                    op = f'INS {dna_sequence}'
                    output_file_sv.write(
                        f'{variant_record.contig} {variant_record.pos} {variant_record.pos} {insertion_prefix}{op}\n')

    output_file_sv.close()
    output_file_snv.close()
    output_file_indel.close()
    print(
        f'Output files generated: {args.output_file_schema}_sv.in, {args.output_file_schema}_snv.in, {args.output_file_schema}_indel.in')
