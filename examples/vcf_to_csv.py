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
    for variant_record in extractor:
        start_chrom = variant_record.contig.replace('chr', '')
        start = variant_record.pos
        ref = variant_record.ref
        alt = variant_record.alt
        length = variant_record.length
        end = variant_record.end
        if variant_record.alt_sv_bracket:
            end_chrom = variant_record.alt_sv_bracket.contig.replace('chr', '')
            if start_chrom != end_chrom:
                end = variant_record.alt_sv_bracket.pos
        else:
            end_chrom = start_chrom

        # Inferred type
        type_inferred = variant_record.variant_type.name
        # Get called EVENTTYPE or SVTYPE from INFO field
        type_called = variant_record.info['EVENTTYPE'] if 'EVENTTYPE' in variant_record.info else \
            variant_record.info['SVTYPE'] if 'SVTYPE' in variant_record.info else None
        # Brackets
        brackets = ''
        if type_inferred == VariantType.DEL.name:
            brackets = 'N['
        elif type_inferred == VariantType.DUP.name:
            brackets = ']N'
        elif type_inferred == VariantType.INV.name:
            prefix = 'N' if variant_record.alt_sv_bracket.prefix else ''
            suffix = 'N' if variant_record.alt_sv_bracket.suffix else ''
            brackets = prefix + variant_record.alt_sv_bracket.bracket + suffix
        elif type_inferred == VariantType.TRN.name:
            prefix = 'N' if variant_record.alt_sv_bracket.prefix else ''
            suffix = 'N' if variant_record.alt_sv_bracket.suffix else ''
            brackets = prefix + variant_record.alt_sv_bracket.bracket + variant_record.alt_sv_bracket.bracket + suffix

        variants.append([start_chrom, start, end_chrom, end, ref, alt, length,
                        brackets, type_inferred, type_called, str(variant_record)])

    df = pd.DataFrame(variants, columns=['start_chrom', 'start', 'end_chrom', 'end', 'ref',
                      'alt', 'length', 'brackets', 'type_inferred', 'type_called', 'variant_record_str'])

    df.drop(['variant_record_str'], axis=1).to_csv(f'{args.output_file}', index=False)
