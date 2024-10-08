# VariantExtractor<!-- omit in toc -->

[![DOI](https://zenodo.org/badge/476634824.svg)](https://zenodo.org/doi/10.5281/zenodo.12755170)

R. Martín et al., “ONCOLINER: A new solution for monitoring, improving, and harmonizing somatic variant calling across genomic oncology centers,” _Cell Genomics_, vol. 4, no. 9. Elsevier BV, p. 100639, Sep. 2024. [doi: 10.1016/j.xgen.2024.100639](https://doi.org/10.1016/j.xgen.2024.100639)

**Deterministic and standard extractor of indels, SNVs and structural variants (SVs)** from VCF files built under the frame of [EUCANCan](https://eucancan.com/)'s second work package. VariantExtractor is a Python package (**requires Python version 3.6 or higher**) and provides a set of data structures and functions to extract variants from VCF files in a **deterministic and standard** way while [adding information](#variantrecord) to facilitate afterwards processing. It homogenizes [multiallelic variants](#multiallelic-variants), [MNPs](#snvs) and [SVs](#structural-variants) (including [imprecise paired breakends](#imprecise-paired-breakends) and [single breakends](#single-breakends)). The package is designed to be used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis. Check the [available documentation](https://eucancan.github.io/variant-extractor/) for more information.

While there is somewhat of an agreement on how to label the SNVs and indels variants, this is not the case for the structural variants. In the current scenario, different labeling between variant callers makes comparisons between structural variants difficult. This package provides an unified interface to extract variants (included structural variants) from VCFs generated by different variant callers. Apart from reading the VCF file, VariantExtractor **adds a preprocessing layer to homogenize the variants** extracted from the file. This way, the variants can be used in downstream analysis in a consistent way. For more information about the homogenization process, check the [homogenization rules](#homogenization-rules) section.


## Table of contents<!-- omit in toc -->
- [Getting started](#getting-started)
  - [Installation](#installation)
- [Usage](#usage)
- [VariantRecord](#variantrecord)
  - [VariantType](#varianttype)
  - [BreakendSVRecord](#breakendsvrecord)
  - [ShorthandSVRecord](#shorthandsvrecord)
- [Homogenization rules](#homogenization-rules)
  - [Multiallelic variants](#multiallelic-variants)
  - [SNVs](#snvs)
  - [Structural variants](#structural-variants)
    - [Breakend vs shorthand notation](#breakend-vs-shorthand-notation)
    - [Paired breakends](#paired-breakends)
    - [Inferred breakend pairs](#inferred-breakend-pairs)
    - [Imprecise paired breakends](#imprecise-paired-breakends)
    - [Single breakends](#single-breakends)
- [Dependencies](#dependencies)


## Getting started
### Installation
VariantExtractor is available on PyPI and can be installed using `pip`:
```bash
pip install variant-extractor
```

## Usage
```python
# Import the package
from variant_extractor import VariantExtractor

# Create a new instance of the class
extractor = VariantExtractor('/path/to/file.vcf')
# Iterate through the variants
for variant_record in extractor:
    print(f'Found variant of type {variant_record.variant_type.name}: {variant_record.contig}:{variant_record.pos}')
```

```python
# Import the package
from variant_extractor import VariantExtractor

# Create a new instance of the class
extractor = VariantExtractor('/path/to/file.vcf')

# Save variants to a CSV file
extractor.to_dataframe().drop(['variant_record_obj'], axis=1).to_csv('/path/to/output.csv', index=False)
```

For a more complete list of examples, check the [examples](./examples/) directory. This folder also includes an example of a [script for normalizing VCF files](examples/normalize_vcf.py) following the [homogenization rules](#homogenization-rules).

## VariantRecord
The `VariantExtractor` constructor returns a generator of `VariantRecord` instances. The `VariantRecord` class is a container for the information contained in a VCF record plus some extra useful information.

| Property           | Type                                                    | Description                                                                                                   |
| ------------------ | ------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| `contig`           | `str`                                                   | Contig name                                                                                                   |
| `pos`              | `int`                                                   | Position on the contig                                                                                        |
| `end`              | `int`                                                   | End position of the variant in the contig (same as `pos` for TRA and SNV)                                     |
| `length`           | `int`                                                   | Length of the variant                                                                                         |
| `id`               | `Optional[str]`                                         | Record identifier                                                                                             |
| `ref`              | `str`                                                   | Reference sequence                                                                                            |
| `alt`              | `str`                                                   | Alternative sequence                                                                                          |
| `qual`             | `Optional[float]`                                       | Quality score for the assertion made in ALT                                                                   |
| `filter`           | `List[str]`                                             | Filter status. `PASS` if this position has passed all filters. Otherwise, it contains the filters that failed |
| `info`             | `Dict[str, Any]`                                        | Additional information                                                                                        |
| `format`           | `List[str]`                                             | Specifies data types and order of the genotype information                                                    |
| `samples`          | `Dict[str, Dict[str, Any]]`                             | Genotype information for each sample                                                                          |
| `variant_type`     | [`VariantType`](#varianttype)                           | Variant type inferred                                                                                         |
| `alt_sv_breakend`  | `Optional[`[`BreakendSVRecord`](#brekendsvrecord)`]`    | Breakend SV info, present only for SVs with breakend notation. For example, `G]17:198982]`                    |
| `alt_sv_shorthand` | `Optional[`[`ShorthandSVRecord`](#shorthandsvrecord)`]` | Shorthand SV info, present only for SVs with shorthand notation. For example, `<DUP:TANDEM>`                  |

### VariantType
The `VariantType` enum describes the type of the variant. For structural variants, it is inferred **only** from the breakend notation (or shorthand notation). It does not take into account any `INFO` field (`SVTYPE` nor `EVENTYPE`) that might be added by the variant caller afterwards.

| REF  | ALT                                      | Variant name | Description                                                          |
| ---- | ---------------------------------------- | ------------ | -------------------------------------------------------------------- |
| A    | G                                        | SNV          | Single nucleotide variant                                            |
| AGTG | A                                        | DEL          | Deletion                                                             |
| A    | A[1:20[ or \<DEL\>                       | DEL          | Deletion                                                             |
| A    | ACCT or \<INS\>                          | INS          | Insertion                                                            |
| A    | ]1:20]A or \<DUP\>                       | DUP          | Duplication                                                          |
| A    | A]1:20] or [1:20[A                       | INV          | Inversion. **[\<INV\> is a special case](#the-special-case-of-inv)** |
| A    | \<CNV\>                                  | CNV          | Copy number variation                                                |
| A    | A]X:20] or A[X:20[ or ]X:20]A or [X:20[A | TRA          | Translocation                                                        |
| A    | A. or .A                                 | SGL          | Single breakend                                                      |

### BreakendSVRecord
The `BreakendSVRecord` class is a container for the information contained in a VCF record for SVs with breakend notation.

| Property  | Type            | Description                                                                                                     |
| --------- | --------------- | --------------------------------------------------------------------------------------------------------------- |
| `prefix`  | `Optional[str]` | Prefix of the SV record with breakend notation. For example, for `G]17:198982]` the prefix will be `G`          |
| `bracket` | `str`           | Bracket of the SV record with breakend notation. For example, for `G]17:198982]` the bracket will be `]`        |
| `contig`  | `str`           | Contig of the SV record with breakend notation. For example, for `G]17:198982]` the contig will be `17`         |
| `pos`     | `int`           | Position of the SV record with breakend notation. For example, for `G]17:198982]` the position will be `198982` |
| `suffix`  | `Optional[str]` | Suffix of the SV record with breakend notation. For example, for `G]17:198982]` the suffix will be `None`       |

### ShorthandSVRecord
The `ShorthandSVRecord` class is a container for the information contained in a VCF record for SVs with shorthand notation.

| Property | Type        | Description                                                                                                                                                                |
| -------- | ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `type`   | `str`       | Type of the SV record with shorthand notation. One of the following, `'DEL'`, `'INS'`, `'DUP'`, `'INV'` or `'CNV'`. For example, for `<DUP:TANDEM>` the type will be `DUP` |
| `extra`  | `List[str]` | Extra information of the SV. For example, for `<DUP:TANDEM:AA>` the extra will be `['TANDEM', 'AA']`                                                                       |

## Homogenization rules
VariantExtractor provides a unified interface to extract variants (included structural variants) from VCF files generated by different variant callers. The variants are homogenized and returned applying the following rules:

### Multiallelic variants
An entry with multiple `ALT` sequences (multiallelic) is divided into multiple entries with a single `ALT` field. This entries with a single `ALT` field are then processed with the rest of the homogeneization rules. For example:

| CHROM | POS | ID             | REF | ALT | FILTER |
| ----- | --- | -------------- | --- | --- | ------ |
| 2     | 1   | multiallelic_1 | G   | C,T | PASS   |

is returned as:

| CHROM | POS | ID               | REF | ALT | FILTER | [`VariantType`](#varianttype) |
| ----- | --- | ---------------- | --- | --- | ------ | ----------------------------- |
| 2     | 1   | multiallelic_1_0 | G   | C   | PASS   | SNV                           |
| 2     | 1   | multiallelic_1_1 | G   | T   | PASS   | SNV                           |


### SNVs
Entries with `REF/ALT` of the same lenghts are treated like SNVs. If the `REF/ALT` sequences are more than one nucleotide (MNPs), they are divided into multiple atomic SNVs. For example:

| CHROM | POS | ID    | REF | ALT | FILTER |
| ----- | --- | ----- | --- | --- | ------ |
| 2     | 1   | snv_1 | C   | G   | PASS   |
| 2     | 3   | mnp_1 | TAG | AGT | PASS   |

are returned as:

| CHROM | POS | ID      | REF | ALT | FILTER | [`VariantType`](#varianttype) |
| ----- | --- | ------- | --- | --- | ------ | ----------------------------- |
| 2     | 3   | snv_1   | C   | G   | PASS   | SNV                           |
| 2     | 3   | mnp_1_0 | T   | A   | PASS   | SNV                           |
| 2     | 4   | mnp_1_1 | A   | G   | PASS   | SNV                           |
| 2     | 5   | mnp_1_2 | G   | T   | PASS   | SNV                           |

<!-- ### Compound indels
All entries with the `REF/ALT` of different lengths are treated as compound indels (or complex indels). They are left-trimmed and divided into multiple atomic SNVs and an insertion (INS) or a deletion (DEL). If the `REF` sequence is longer than the `ALT` sequence, it is considered a deletion. If the `REF` sequence is shorter than the `ALT` sequence, it is considered an insertion. For example:

| CHROM | POS  | ID           | REF     | ALT       | FILTER |
| ----- | ---- | ------------ | ------- | --------- | ------ |
| 1     | 2000 | standard_del | CT      | C         | PASS   |
| 1     | 2100 | standard_ins | C       | CAA       | PASS   |
| 1     | 2200 | compund_del  | CCTGAAA | CGA       | PASS   |
| 1     | 2300 | compund_ins  | GT      | CAATATATA | PASS   |

are returned as:

| CHROM | POS  | ID            | REF   | ALT      | FILTER | [`VariantType`](#varianttype) |
| ----- | ---- | ------------- | ----- | -------- | ------ | ----------------------------- |
| 1     | 2000 | standard_del  | CT    | C        | PASS   | DEL                           |
| 1     | 2100 | standard_ins  | C     | CAA      | PASS   | INS                           |
| 1     | 2201 | compund_del_0 | C     | G        | PASS   | SNV                           |
| 1     | 2202 | compund_del_1 | TGAAA | T        | PASS   | DEL                           |
| 1     | 2202 | compund_del_2 | T     | A        | PASS   | SNV                           |
| 1     | 2300 | compund_ins_0 | G     | C        | PASS   | SNV                           |
| 1     | 2301 | compund_ins_1 | T     | TATATATA | PASS   | DEL                           |
| 1     | 2301 | compund_ins_2 | T     | A        | PASS   | SNV                           | --> |


### Structural variants
VariantExtractor returns one entry per structural variant (one entry per breakend pair). This helps to avoid the ambiguity of the notation and keeps the process deterministic. For this reason, in case of paired breakends, the breakend with the lowest chromosome and/or position is returned. If a breakend is not the lowest chromosome and/or position and is missing its pair, its pair is [inferred and returned](#inferred-breakend-pairs).

#### Breakend vs shorthand notation
Entries with the same information, either described with shorthand or breakend notation, will be returned the same way. Here is an example for a DEL entry:

| CHROM | POS  | ID        | REF         | ALT       | FILTER | INFO                 |
| ----- | ---- | --------- | ----------- | --------- | ------ | -------------------- |
| 1     | 3000 | event_1_o | A           | A[1:5000[ | PASS   | SVTYPE=BND           |
| 1     | 5000 | event_1_h | A           | ]1:3000]A | PASS   | SVTYPE=BND           |
| 1     | 3000 | event_1   | A           | A[1:5000[ | PASS   | SVTYPE=DEL           |
| 1     | 3000 | event_1   | A           | \<DEL\>   | PASS   | SVTYPE=DEL; END=5000 |
| 1     | 3000 | event_1   | AGTCACAA... | A         | PASS   |                      |

are returned as one entry (each one of them with their own `ALT` field), but with the same `VariantRecord.end` and `VariantType`:

| CHROM | POS  | ID      | REF | ALT | FILTER | INFO | [`VariantType`](#varianttype) | [`VariantRecord.end`](#variantrecord) |
| ----- | ---- | ------- | --- | --- | ------ | ---- | ----------------------------- | ------------------------------------- |
| 1     | 3000 | event_1 | A   | ... | PASS   | ...  | DEL                           | 5000                                  |

##### The special case of INV<!-- omit in toc -->
\<INV\> is a special case of shorthand notation because it represents two paired breakends. For example, the following shorthand notation:

| CHROM | POS    | ID      | REF | ALT     | FILTER | INFO                  |
| ----- | ------ | ------- | --- | ------- | ------ | --------------------- |
| 2     | 321682 | event_1 | T   | \<INV\> | PASS   | SVTYPE=INV;END=421681 |

is equivalent to the following breakends:

| CHROM | POS    | ID        | REF | ALT         | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ------ | --------- | --- | ----------- | ------ | ---------- | ----------------------------- |
| 2     | 321681 | event_1_0 | N   | N]2:421681] | PASS   | SVTYPE=INV | INV                           |
| 2     | 321682 | event_1_1 | T   | [2:421682[T | PASS   | SVTYPE=INV | INV                           |

In this case, VariantExtractor converts internally \<INV\> to two entries with breakend notation (one for each breakend pair). Note that the `N` will be replaced with the correct nucleotide if `fasta_ref` is provided to VariantExtractor.


#### Paired breakends
For **paired breakends**, breakends are paired using the `INFO` fields `MATEID` or `PARID`. If these fields are not available, they are paired using their coordinates (contig+position). The breakend with the lowest chromosome and/or position is returned. For example:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- |
| 2     | 3000 | event_1_o | T   | ]3:5000]T | PASS   | SVTYPE=BND |
| 3     | 5000 | event_1_h | G   | G[2:3000[ | PASS   | SVTYPE=BND |
| 1     | 3000 | event_2_o | A   | A[1:5000[ | PASS   | SVTYPE=BND |
| 1     | 5000 | event_2_h | A   | ]1:3000]A | PASS   | SVTYPE=BND |

are returned as one entry per variant:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- | ----------------------------- |
| 2     | 3000 | event_1_o | T   | ]3:5000]T | PASS   | SVTYPE=BND | TRA                           |
| 1     | 3000 | event_2_o | A   | A[1:5000[ | PASS   | SVTYPE=BND | DEL                           |


#### Inferred breakend pairs
If **all** breakends are missing their pair, the breakends with the lowest chromosome and/or position are inferred and returned. For example:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- |
| 3     | 5000 | event_1_h | G   | G[2:3000[ | PASS   | SVTYPE=BND |
| 1     | 5000 | event_2_h | A   | ]1:3000]A | PASS   | SVTYPE=BND |

are returned as their inferred breakend pair with the lowest chromosome and/or position:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- | ----------------------------- |
| 2     | 3000 | event_1_h | N   | ]3:5000]N | PASS   | SVTYPE=BND | TRA                           |
| 1     | 3000 | event_2_h | A   | A[1:5000[ | PASS   | SVTYPE=BND | DEL                           |

Note that the `N` will be replaced with the correct nucleotide if `fasta_ref` is provided to VariantExtractor. The following equivalencies are applied:

| CHROM1 | POS1 | REF1 | ALT1     | CHROM2 | POS2 | REF2 | ALT2     |
| ------ | ---- | ---- | -------- | ------ | ---- | ---- | -------- |
| 1      | 500  | N    | N[7:800[ | 7      | 800  | N    | ]1:500]N |
| 1      | 500  | N    | ]7:800]N | 7      | 800  | N    | N[1:500[ |
| 1      | 500  | N    | [7:800[N | 7      | 800  | N    | [1:500[N |
| 1      | 500  | N    | N]7:800] | 7      | 800  | N    | N]1:500] |
      

#### Imprecise paired breakends
Imprecise breakends do not match exactly with their pair in coordinates. In this case, they are paired using the `INFO` fields `MATEID` or `PARID`. As with the rest of variants, for each breakend pair, only the breakend with the lowest chromosome and/or position is returned. However, it is important to notice that the `CIPOS` field is lost for the other breakend. For example:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO                                   |
| ----- | ---- | --------- | --- | --------- | ------ | -------------------------------------- |
| 2     | 3010 | event_1_o | T   | T[3:5000[ | PASS   | SVTYPE=BND;CIPOS=0,50;PARID=event_1_h  |
| 3     | 5050 | event_1_h | A   | ]2:3050]A | PASS   | SVTYPE=BND;CIPOS=0,100;PARID=event_1_o |

are paired and the entry with the lowest chromosome and/or position is returned:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO                            | [`VariantType`](#varianttype) |
| ----- | ---- | --------- | --- | --------- | ------ | ------------------------------- | ----------------------------- |
| 2     | 3010 | event_1_o | T   | T[3:5000[ | PASS   | SVTYPE=BND;CIPOS=0,50;PARID=a_h | TRA                           |


#### Single breakends
Single breakends cannot be matched with other breakends because they lack a mate. They may be able to be matched later in downstream analysis. That is why each one is kept as a different variant. For example:

| CHROM | POS  | ID      | REF | ALT | FILTER | INFO       |
| ----- | ---- | ------- | --- | --- | ------ | ---------- |
| 2     | 3000 | event_s | T   | T.  | PASS   | SVTYPE=BND |
| 3     | 5000 | event_m | G   | .G  | PASS   | SVTYPE=BND |

are returned as two entries:

| CHROM | POS  | ID      | REF | ALT | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ---- | ------- | --- | --- | ------ | ---------- | ----------------------------- |
| 2     | 3000 | event_s | T   | T.  | PASS   | SVTYPE=BND | SGL                           |
| 3     | 5000 | event_m | G   | .G  | PASS   | SVTYPE=BND | SGL                           |


## Dependencies

The dependencies are covered by their own respective licenses as follows:

* [Python/Pysam package](https://github.com/pysam-developers/pysam) (MIT license)
