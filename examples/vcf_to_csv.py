# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
'''
Generates a CSV file from an input VCF file
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
    from variant_extractor import VariantExtractor
    from variant_extractor.variants import VariantType

    # Parse arguments
    parser = ArgumentParser(description='Generate CSV file from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_file', help='Output file')
    parser.add_argument('-f', '--fasta-ref', help='FASTA reference file')
    args = parser.parse_args()

    variants = []

    print(f'Reading VCF file: {args.vcf_file}')
    extractor = VariantExtractor(args.vcf_file, fasta_ref=args.fasta_ref)
    extractor.to_dataframe().drop(['variant_record_obj'], axis=1).to_csv(f'{args.output_file}', index=False)
