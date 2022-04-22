# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
'''
Generates a normalized VCF file from a VCF
Expected usage:
    $ python normalize.py <vcf_file> <output_vcf_file>
Use --help for more information.
'''
from argparse import ArgumentParser


def _extract_header(vcf_file):
    header = ''
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header += line
            else:
                break
    return header


if __name__ == '__main__':
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + '/../src/')
    from variant_extractor import VariantExtractor
    from variant_extractor.variants import VariantType

    # Parse arguments
    parser = ArgumentParser(description='Generate normalized VCF file from a VCF file')
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('output_vcf_file', help='Output VCF file')
    args = parser.parse_args()

    if not args.vcf_file.endswith('.vcf'):
        raise ValueError('Input file must be a VCF file')

    # Open output file, write as stream
    with open(args.output_vcf_file, 'w') as output_vcf:
        # Write header
        output_vcf.write(_extract_header(args.vcf_file))
        print(f'Reading {args.vcf_file}...')
        # Open input file, read with variant_extractor
        extractor = VariantExtractor()
        for variant_record in extractor.read_vcf(args.vcf_file):
            output_vcf.write(str(variant_record)+'\n')
