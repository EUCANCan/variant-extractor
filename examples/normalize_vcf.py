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


def _convert_info_key_value(key, value):
    if value is None:
        return key
    elif isinstance(value, str):
        return f'{key}={value}'
    elif hasattr(value, '__len__'):
        return key+'=' + ','.join([str(v) for v in value])
    else:
        return key+'='+str(value)


def _convert_sample_value(key, value):
    if key == 'GT':
        return '/'.join([str(v) if v is not None else '.' for v in value])
    elif value is None:
        return '.'
    elif isinstance(value, str):
        return value
    elif hasattr(value, '__len__'):
        return ','.join([str(v) for v in value])
    else:
        return str(value)


def _convert_to_vcf(variant_record):
    contig = variant_record.contig
    pos = variant_record.pos
    id_ = variant_record.id if variant_record.id else '.'
    ref = variant_record.ref
    alt = variant_record.alt
    qual = variant_record.qual if variant_record.qual else '.'
    filter_ = ";".join(variant_record.filter) if variant_record.filter else '.'
    info = ";".join([_convert_info_key_value(k, v) for k, v in variant_record.info.items()])
    format_ = ":".join(variant_record.format)
    samples_list = [":".join([_convert_sample_value(k, v) for k, v in variant_record.samples[sample_name].items()])
                    for sample_name in variant_record.samples]
    samples = "\t".join(samples_list)
    return f'{contig}\t{pos}\t{id_}\t{ref}\t{alt}\t{qual}\t{filter_}\t{info}\t{format_}\t{samples}\n'


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
            output_vcf.write(_convert_to_vcf(variant_record))
