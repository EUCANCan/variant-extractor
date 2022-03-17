# Variant extractor from VCF files
Extractor of INDELs, SNVs and structural variations (SVs) from VCF files built under the frame of [EUCANCan](https://eucancan.com/)'s second work package. This package provides a Python (requires Python version 3.6 or higher) package called `variant_extractor` with a set of classes and functions to extract variants from VCF files.

DNA variants refer to differences among the genomes of members of the same species and consist in single (SNV), small (INDELs) or large structural differences (SVs). However, while there is somewhat of an agreement on how to label the SNVs and INDELs variants, this is not the case for the structural variants. In the current scenario of structural variants, different labeling between variant callers makes comparisons between results difficult. This one of the reasons why there is a need to standardize variant labeling procedures for cancer genomic analysis. This is precisely one of the objectives of EUCANCan's second work package.


## Features
The objective of this **variant extractor** is to ingest different VCF files from different variant callers and extract the variants in a standardized and deterministic format.

* Returns only one entry per variant. For paired breakends, since they are redudant, the breakend with the lowest chromosome and/or position is kept. For example:

| CHROM | POS  | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --- | --------- | ------ | ---------- |
| 2     | 3000 | T   | T[3:5000[ | PASS   | SVTYPE=BND |
| 3     | 5000 | T   | ]2:3000]T | PASS   | SVTYPE=BND |
| 1     | 3000 | A   | A[1:5000[ | PASS   | SVTYPE=BND |
| 1     | 5000 | A   | ]1:3000]A | PASS   | SVTYPE=BND |

Would be treated as:

| CHROM | POS  | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --- | --------- | ------ | ---------- |
| 2     | 3000 | T   | T[3:5000[ | PASS   | SVTYPE=TRN |
| 1     | 3000 | T   | T[1:5000[ | PASS   | SVTYPE=DEL |


* The variation type for structural variations is inferred **only** from the bracket notation, it does not take into account the `INFO` fields `"SVTYPE"` or `"EVENTYPE"` (this is pending manual validation).

* Entries with `REF/ALT` of the same lenghts are treated like multiple SNVs. For example:

| CHROM | POS | REF | ALT | FILTER |
| ----- | --- | --- | --- | ------ |
| 2     | 3   | TAG | AGT | PASS   |

Would be treated as:

| CHROM | POS | REF | ALT | FILTER |
| ----- | --- | --- | --- | ------ |
| 2     | 3   | T   | A   | PASS   |
| 2     | 4   | A   | G   | PASS   |
| 2     | 5   | G   | T   | PASS   |

## Caveats

* Right now, VCF entries with multiple ALT values are ignored.

* Imprecise structural variants breakends with bracket notation are paired using the `INFO` fields `"MATEID"` or `"PARID"`. As with the rest of variants, only the breakend with the lowest chromosome and/or position is kept. This way information is lost about the uncertainty (`"CIPOS"`) in the other breakend, but it keeps a one entry per variant standard.

## Getting started

### Examples
