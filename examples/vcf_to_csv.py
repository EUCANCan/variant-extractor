# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# MIT License
'''
Generates a .CSV input from a VCF file
Expected usage:
    $ python vcf_to_csv.py <vcf_file> <output_file>
Use --help for more information.
'''
import pandas as pd
from argparse import ArgumentParser

if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + '/../src/')
    from variant_extractor import VariantExtractor, VariationType  

    # Parse arguments
    parser = ArgumentParser(description='Generate CSV file from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_file', help='Output file')
    args = parser.parse_args()

    variants = []

    print(f'Reading VCF file: {args.vcf_file}')
    extractor = VariantExtractor()
    for var_type, variant_record in extractor.read_vcf(args.vcf_file):
        start_chrom = variant_record.contig.replace('chr', '')
        start = variant_record.pos
        ref = variant_record.ref
        alt = variant_record.alts[0]
        length = variant_record.end - variant_record.pos
        end = variant_record.end
        if variant_record.alt_sv_bracket:
            end_chrom = variant_record.alt_sv_bracket.contig.replace('chr', '')
            if start_chrom != end_chrom:
                end = variant_record.alt_sv_bracket.pos
        else:
            end_chrom = variant_record.contig.replace('chr', '')
        # Hotfix for indels INS
        if not variant_record.alt_sv_bracket and not variant_record.alt_sv_shorthand and var_type == VariationType.INS or var_type == VariationType.INDEL_INS:
            length = len(alt)-len(ref)
        variants.append([start_chrom, start, end_chrom, end, ref, alt, length, var_type.name])

    df = pd.DataFrame(variants, columns=['start_chrom', 'start',
                      'end_chrom', 'end', 'ref', 'alt', 'length', 'var_type'])
    df.to_csv(f'{args.output_file}.csv', index=False)
