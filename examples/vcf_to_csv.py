# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
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
    from variant_extractor import VariantExtractor
    from variant_extractor.variants import VariantType

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
            if var_type == VariantType.INS and 'SVLEN' in variant_record.info:
                if isinstance(variant_record.info['SVLEN'], tuple):
                    length = int(variant_record.info['SVLEN'][0])
                else:
                    length = int(variant_record.info['SVLEN'])
        # Hotfix for indels INS
        if not variant_record.alt_sv_bracket and not variant_record.alt_sv_shorthand and var_type == VariantType.INS:
            length = len(alt)-len(ref)

        # Inferred type
        type_inferred = var_type.name
        # Get called EVENTTYPE or SVTYPE from INFO field
        type_called = variant_record.info['EVENTTYPE'] if 'EVENTTYPE' in variant_record.info else \
            variant_record.info['SVTYPE'] if 'SVTYPE' in variant_record.info else None
        variants.append([start_chrom, start, end_chrom, end, ref, alt, length, type_inferred, type_called])

    df = pd.DataFrame(variants, columns=['start_chrom', 'start',
                      'end_chrom', 'end', 'ref', 'alt', 'length', 'type_inferred', 'type_called'])
    df.to_csv(f'{args.output_file}', index=False)
