# VCF Variant extractor <!-- omit in toc -->
**Extractor of INDELs, SNVs and structural variants (SVs)** from VCF files built under the frame of [EUCANCan](https://eucancan.com/)'s second work package. While there is somewhat of an agreement on how to label the SNVs and INDELs variants, this is not the case for the structural variants. In the current scenario of structural variants, different labeling between variant callers makes comparisons between results difficult. This one of the reasons why there is a need to standardize variant labeling procedures for cancer genomic analysis. This is precisely one of the objectives of EUCANCan's second work package.

This package provides a Python (**requires Python version 3.6 or higher**) package called `variant_extractor` with a set of classes and functions to read and extract variants from VCF files. The package is designed to be used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis.

## Table of contents <!-- omit in toc -->
- [Installation](#installation)
- [Getting started](#getting-started)
  - [Examples](#examples)
- [How it works](#how-it-works)
  - [SNVs](#snvs)
  - [Indels](#indels)
  - [Structural variants](#structural-variants)
- [Caveats](#caveats)

## Installation
(Not yet uploaded to pypi.org)
```
pip install variant_extractor
```

## Getting started
```python
# Import the package
from variant_extractor import VariantExtractor, variantType  

# Create a new instance of the class
extractor = VariantExtractor()
# Read the VCF file
variants = extractor.read_vcf('/path/to/file.vcf')
# Extract the variants
for variant_type, variant_record in variants:
    print(f'Found variant of type {variant_type.name()}: {variant_record.contig}:{variant_record.pos}')
```

### Examples
For a list of examples, see the [examples](./examples/) directory.


## How it works
The objective of this **variant extractor** is to ingest different VCF files from different variant callers and extract the variants in a standardized and deterministic format. This way, it is easier to compare the results of different variant callers.

In order to achieve this goal, the following steps are performed for each variant type:
### SNVs
Entries with `REF/ALT` of the same lenghts are treated like multiple SNVs. For example:

| CHROM | POS | REF | ALT | FILTER |
| ----- | --- | --- | --- | ------ |
| 2     | 3   | TAG | AGT | PASS   |

Would be treated as:

| CHROM | POS | REF | ALT | FILTER | VAR_TYPE |
| ----- | --- | --- | --- | ------ | -------- |
| 2     | 3   | T   | A   | PASS   | SNV      |
| 2     | 4   | A   | G   | PASS   | SNV      |
| 2     | 5   | G   | T   | PASS   | SNV      |

### Indels
All entries with the following format would be treated as indels:

| CHROM | POS | REF | ALT | FILTER | VAR_TYPE  |
| ----- | --- | --- | --- | ------ | --------- |
| 2     | 3   | T   | TGT | PASS   | INDEL_INS |
| 2     | 10  | TA  | AAT | PASS   | INDEL_INS |
| 2     | 5   | TAG | T   | PASS   | INDEL_DEL |
| 2     | 8   | TAG | AG  | PASS   | INDEL_DEL |

### Structural variants
The variant type for structural variants is inferred **only** from the bracket notation, it does not take into account the `INFO` fields `"SVTYPE"` or `"EVENTYPE"` (this is pending manual validation).

Returns only one entry per variant. For **paired breakends**, since they are redudant, the breakend with the lowest chromosome and/or position is kept. For example:

| CHROM | POS  | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --- | --------- | ------ | ---------- |
| 2     | 3000 | T   | T[3:5000[ | PASS   | SVTYPE=BND |
| 3     | 5000 | T   | ]2:3000]T | PASS   | SVTYPE=BND |
| 1     | 3000 | A   | A[1:5000[ | PASS   | SVTYPE=BND |
| 1     | 5000 | A   | ]1:3000]A | PASS   | SVTYPE=BND |

Would be treated as:

| CHROM | POS  | REF | ALT       | FILTER | INFO       | VAR_TYPE |
| ----- | ---- | --- | --------- | ------ | ---------- | -------- |
| 2     | 3000 | T   | T[3:5000[ | PASS   | SVTYPE=BND | TRN      |
| 1     | 3000 | A   | A[1:5000[ | PASS   | SVTYPE=BND | DEL      |

Returns an entry with the same information for both `<>` notation and bracket notation. For example:

| CHROM | POS  | REF | ALT       | FILTER | INFO                 |
| ----- | ---- | --- | --------- | ------ | -------------------- |
| 1     | 3000 | A   | A[1:5000[ | PASS   | SVTYPE=BND           |
| 1     | 5000 | A   | ]1:3000]A | PASS   | SVTYPE=BND           |
| 1     | 3000 | A   | A[1:5000[ | PASS   | SVTYPE=DEL           |
| 1     | 3000 | A   | `<DEL>`   | PASS   | SVTYPE=DEL; END=5000 |

Would be all treated as the same entry:

| CHROM | POS  | REF | ALT | FILTER | INFO | VAR_TYPE | END  |
| ----- | ---- | --- | --- | ------ | ---- | -------- | ---- |
| 2     | 3000 | A   | ... | PASS   | ...  | DEL      | 5000 |


## Caveats

* VCF entries with multiple ALT values are ignored.

* Imprecise structural variants breakends with bracket notation are paired using the `INFO` fields `"MATEID"` or `"PARID"`. As with the rest of variants, only the breakend with the lowest chromosome and/or position is kept. This way information is lost about the uncertainty (`"CIPOS"`) in the other breakend, but it keeps a one entry per variant standard.
