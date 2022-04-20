# VCF variant extractor<!-- omit in toc -->
**Deterministic and standard extractor of indels, SNVs and structural variants (SVs)** from VCF files built under the frame of [EUCANCan](https://eucancan.com/)'s second work package. `variant_extractor` is a Python package (**requires Python version 3.6 or higher**) and provides a set of data structures and functions to extract variants from VCF files in a **deterministic and standard** way while [adding information](#variantrecord) to facilitate afterwards processing. It homogenizes [multiallelic variants](#multiallelic-variants), [MNPs](#snvs), [compound indels](#compound-indels) and [SVs](#svs) (including [imprecise paired breakends](#imprecise-paired-breakends) and [single breakends](#single-breakends)). The package is designed to be used in a pipeline, where the variants are ingested from VCF files and then used in downstream analysis. Check the [available documentation](https://eucancan.github.io/variant-extractor/) for more information.

While there is somewhat of an agreement on how to label the SNVs and indels variants, this is not the case for the structural variants. In the current scenario, different labeling between variant callers makes comparisons between structural variants difficult. This package provides an unified interface to extract variants (included structural variants) from VCFs generated by different variant callers. Apart from reading the VCF file, the `variant_extractor` **adds a preprocessing layer to homogenize the variants** extracted from the file. This way, the variants can be used in downstream analysis in a consistent way. For more information about the homogenization process, check the [homogenization rules](#homogenization-rules) section.


## Table of contents<!-- omit in toc -->
- [Usage](#usage)
- [VariantRecord](#variantrecord)
  - [BracketSVRecord](#bracketsvrecord)
  - [ShorthandSVRecord](#shorthandsvrecord)
- [VariantType](#varianttype)
- [Homogenization rules](#homogenization-rules)
  - [Multiallelic variants](#multiallelic-variants)
  - [SNVs](#snvs)
  - [Compound indels](#compound-indels)
  - [Structural variants](#structural-variants)
    - [Bracket vs shorthand notation](#bracket-vs-shorthand-notation)
    - [Paired breakends](#paired-breakends)
    - [Inferred pairs](#inferred-pairs)
    - [Imprecise paired breakends](#imprecise-paired-breakends)
    - [Single breakends](#single-breakends)


## Usage
```python
# Import the package
from variant_extractor import VariantExtractor

# Create a new instance of the class
extractor = VariantExtractor()
# Read the VCF file
variants = extractor.read_vcf('/path/to/file.vcf')
# Iterate through the variants
for variant_record in variants:
    print(f'Found variant of type {variant_record.variant_type.name}: {variant_record.contig}:{variant_record.pos}')
```

For a more complete list of examples, check the [examples](./examples/) directory.

## VariantRecord
The `VariantExtractor.read_vcf()` method returns a list of `VariantRecord`. The `VariantRecord` class is a container for the information contained in a VCF record.

| Property           | Type                          | Description                                                                                  |
| ------------------ | ----------------------------- | -------------------------------------------------------------------------------------------- |
| `contig`           | `str`                         | Contig name                                                                                  |
| `pos`              | `int`                         | Position on the contig                                                                       |
| `end`              | `int`                         | End position of the variant in the contig (same as `pos` for TRN and SNV)                    |
| `length`           | `int`                         | Length of the variant                                                                        |
| `id`               | `str`                         | Record identifier                                                                            |
| `ref`              | `str`                         | Reference sequence                                                                           |
| `alt`              | `str`                         | Alternative sequence                                                                         |
| `qual`             | `Optional[float]`             | Quality score for the assertion made in ALT                                                  |
| `filter`           | `pysam.VariantRecordFilter`   | Filter                                                                                       |
| `info`             | `pysam.VariantRecordInfo`     | Dictionary of information fields                                                             |
| `variant_type`     | [`VariantType`](#varianttype) | Variant type inferred                                                                        |
| `alt_sv_bracket`   | `Optional[BracketSVRecord]`   | Bracketed SV info, present only for SVs with bracket notation. For example, `G]17:198982]`   |
| `alt_sv_shorthand` | `Optional[ShorthandSVRecord]` | Shorthand SV info, present only for SVs with shorthand notation. For example, `<DUP:TANDEM>` |

### BracketSVRecord
The `BracketSVRecord` class is a container for the information contained in a VCF record for SVs with bracket notation.

| Property  | Type            | Description                                                                                                    |
| --------- | --------------- | -------------------------------------------------------------------------------------------------------------- |
| `prefix`  | `Optional[str]` | Prefix of the SV record with bracket notation. For example, for `G]17:198982]` the prefix will be `G`          |
| `bracket` | `str`           | Bracket of the SV record with bracket notation. For example, for `G]17:198982]` the bracket will be `]`        |
| `contig`  | `str`           | Contig of the SV record with bracket notation. For example, for `G]17:198982]` the contig will be `17`         |
| `pos`     | `int`           | Position of the SV record with bracket notation. For example, for `G]17:198982]` the position will be `198982` |
| `suffix`  | `Optional[str]` | Suffix of the SV record with bracket notation. For example, for `G]17:198982]` the suffix will be `None`       |

### ShorthandSVRecord
The `ShorthandSVRecord` class is a container for the information contained in a VCF record for SVs with shorthand notation.

| Property | Type        | Description                                                                                                                                                                |
| -------- | ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `type`   | `str`       | Type of the SV record with shorthand notation. One of the following, `'DEL'`, `'INS'`, `'DUP'`, `'INV'` or `'CNV'`. For example, for `<DUP:TANDEM>` the type will be `DUP` |
| `extra`  | `List[str]` | Extra information of the SV. For example, for `<DUP:TANDEM:AA>` the extra will be `['TANDEM', 'AA']`                                                                       |

## VariantType
The `VariantType` enum is a container for the information about the type of the variant. For structural variants is inferred **only** from the bracket notation, it does not take into account any `INFO` (fields `SVTYPE` or `EVENTYPE`).

| REF       | ALT                                      | Variant name | Description                                                   |
| --------- | ---------------------------------------- | ------------ | ------------------------------------------------------------- |
| A         | G                                        | SNV          | Single nucleotide variant                                     |
| AGTG or A | A or A[1:20[ or \<DEL\>                  | DEL          | Deletion                                                      |
| A         | ACCT or \<INS\>                          | INS          | Insertion                                                     |
| A         | ]1:20]a or \<DUP\>                       | DUP          | Duplication                                                   |
| A         | A]1:20] or [1:20[A                       | INV          | Inversion. **[\<INV\> is a special case](#inv-special-case)** |
| A         | \<CNV\>                                  | CNV          | Copy number variation                                         |
| A         | A]X:20] or A[X:20[ or ]X:20]A or [X:20[A | TRN          | Translocation                                                 |
| A         | A. or .A                                 | SGL          | Single breakend                                               |

## Homogenization rules
The `variant_extractor` package provides a unified interface to extract variants (included structural variants) from VCF files generated by different variant callers. The variants are homogenized and returned applying the following rules:

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

| CHROM | POS | ID      | REF | ALT | FILTER |
| ----- | --- | ------- | --- | --- | ------ |
| 2     | 1   | event_1 | C   | G   | PASS   |
| 2     | 3   | event_2 | TAG | AGT | PASS   |

are returned as:

| CHROM | POS | ID        | REF | ALT | FILTER | [`VariantType`](#varianttype) |
| ----- | --- | --------- | --- | --- | ------ | ----------------------------- |
| 2     | 3   | event_1   | C   | G   | PASS   | SNV                           |
| 2     | 3   | event_2_0 | T   | A   | PASS   | SNV                           |
| 2     | 4   | event_2_1 | A   | G   | PASS   | SNV                           |
| 2     | 5   | event_2_2 | G   | T   | PASS   | SNV                           |

### Compound indels
All entries with the `REF/ALT` of different lengths are treated as compound indels. They are left-trimmed and divided into multiple atomic SNVs and an insertion (INS) or a deletion (DEL). If the `REF` sequence is longer than the `ALT` sequence, it is considered a deletion. If the `REF` sequence is shorter than the `ALT` sequence, it is considered an insertion. For example:

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
| 1     | 2301 | compund_ins_2 | T     | A        | PASS   | SNV                           |


### Structural variants
`variant_extractor` returns one entry per structural variant (one entry per breakend pair). This helps to avoid the ambiguity of the notation and keeps the process deterministic. Also, for this reason, in case of paired breakends, the breakend with the lowest chromosome and/or position is returned.

#### Bracket vs shorthand notation
Entries with the same information, either described with shorthand or bracket notation, will be tagged the same way. Here is an example for a DEL entry:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO                 |
| ----- | ---- | --------- | --- | --------- | ------ | -------------------- |
| 1     | 3000 | event_1_o | A   | A[1:5000[ | PASS   | SVTYPE=BND           |
| 1     | 5000 | event_1_h | A   | ]1:3000]A | PASS   | SVTYPE=BND           |
| 1     | 3000 | event_1   | A   | A[1:5000[ | PASS   | SVTYPE=DEL           |
| 1     | 3000 | event_1   | A   | \<DEL\>   | PASS   | SVTYPE=DEL; END=5000 |

are returned as one entry (each one of them with their own `ALT` field), but with the same `VariantRecord.end` and `VariantType`:

| CHROM | POS  | ID      | REF | ALT | FILTER | INFO | [`VariantType`](#varianttype) | [`VariantRecord.end`](#variantrecord) |
| ----- | ---- | ------- | --- | --- | ------ | ---- | ----------------------------- | ------------------------------------- |
| 2     | 3000 | event_1 | A   | ... | PASS   | ...  | DEL                           | 5000                                  |

##### INV special case<!-- omit in toc -->
\<INV\> is a special case of shorthand notation because it represents two paired breakends. For example, the following shorthand notation:

| CHROM | POS    | ID   | REF | ALT     | FILTER | INFO                  |
| ----- | ------ | ---- | --- | ------- | ------ | --------------------- |
| 2     | 321682 | INV0 | T   | \<INV\> | PASS   | SVTYPE=INV;END=421681 |

is equivalent to the following breakends:

| CHROM | POS    | ID      | REF | ALT         | FILTER | INFO       |
| ----- | ------ | ------- | --- | ----------- | ------ | ---------- |
| 2     | 321681 | event_1 | .   | .]2:421681] | PASS   | SVTYPE=BND |
| 2     | 321682 | event_1 | T   | [2:421682[T | PASS   | SVTYPE=BND |

In this case, `variant_extractor` converts internally \<INV\> to two entries with bracket notation (one for each breakend pair).


#### Paired breakends
For **paired breakends**, breakends are paired using their coordinates (contig+position). The breakend with the lowest chromosome and/or position is returned. For example:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- |
| 2     | 3000 | event_1_o | T   | ]3:5000]T | PASS   | SVTYPE=BND |
| 3     | 5000 | event_1_h | G   | G[2:3000[ | PASS   | SVTYPE=BND |
| 1     | 3000 | event_2_o | A   | A[1:5000[ | PASS   | SVTYPE=BND |
| 1     | 5000 | event_2_h | A   | ]1:3000]A | PASS   | SVTYPE=BND |

are returned as one entry per variant:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- | ----------------------------- |
| 2     | 3000 | event_1_o | T   | ]3:5000]T | PASS   | SVTYPE=BND | TRN                           |
| 1     | 3000 | event_2_o | A   | A[1:5000[ | PASS   | SVTYPE=BND | DEL                           |


#### Inferred pairs
If **all** the breakends are missing their pair, the breakend with the lowest chromosome and/or position is inferred and returned. For example:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- |
| 3     | 5000 | event_1_h | G   | G[2:3000[ | PASS   | SVTYPE=BND |
| 1     | 5000 | event_2_h | A   | ]1:3000]A | PASS   | SVTYPE=BND |

are returned as their inferred pair with the lowest chromosome and/or position:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ---- | --------- | --- | --------- | ------ | ---------- | ----------------------------- |
| 2     | 3000 | event_1_o | .   | ]3:5000]. | PASS   | SVTYPE=BND | TRN                           |
| 1     | 3000 | event_2_h | A   | A[1:5000[ | PASS   | SVTYPE=BND | DEL                           |

The following equivalencies are applied:

| CHROM1 | POS1 | REF1 | ALT1     | CHROM2 | POS2 | REF2 | ALT2     |
| ------ | ---- | ---- | -------- | ------ | ---- | ---- | -------- |
| 1      | 500  | N    | N[7:800[ | 7      | 800  | N    | ]1:500]N |
| 1      | 500  | N    | ]7:800]N | 7      | 800  | N    | N[1:500[ |
| 1      | 500  | N    | [7:800[N | 7      | 800  | N    | [1:500[N |
| 1      | 500  | N    | N]7:800] | 7      | 800  | N    | N]1:500] |
      

#### Imprecise paired breakends
Imprecise breakends with bracket notation are paired using the `INFO` fields `MATEID` or `PARID` instead of their coordinates (since they may not match). In order to keep the deterministic process, as with the rest of variants, only the breakend with the lowest chromosome and/or position is returned. However, it is important to notice that the uncertainty information (`CIPOS` field) is lost for the other breakend. For example:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO                             |
| ----- | ---- | --------- | --- | --------- | ------ | -------------------------------- |
| 2     | 3010 | event_1_o | T   | T[3:5000[ | PASS   | SVTYPE=BND;CIPOS=0,50;PARID=a_h  |
| 3     | 5050 | event_1_h | T   | ]2:3050]T | PASS   | SVTYPE=BND;CIPOS=0,100;PARID=a_o |

are returned as one entry:

| CHROM | POS  | ID        | REF | ALT       | FILTER | INFO                            | [`VariantType`](#varianttype) |
| ----- | ---- | --------- | --- | --------- | ------ | ------------------------------- | ----------------------------- |
| 2     | 3010 | event_1_o | T   | T[3:5000[ | PASS   | SVTYPE=BND;CIPOS=0,50;PARID=a_h | TRN                           |


#### Single breakends
Single breakends cannot be matched with other breakends because they lack a mate. That is why each one is kept as a different variant. For example:

| CHROM | POS  | ID  | REF | ALT | FILTER | INFO       |
| ----- | ---- | --- | --- | --- | ------ | ---------- |
| 2     | 3000 | a_s | T   | T.  | PASS   | SVTYPE=BND |
| 3     | 5000 | m_s | G   | .G  | PASS   | SVTYPE=BND |

are be returned as two entries:

| CHROM | POS  | ID  | REF | ALT | FILTER | INFO       | [`VariantType`](#varianttype) |
| ----- | ---- | --- | --- | --- | ------ | ---------- | ----------------------------- |
| 2     | 3000 | a_s | T   | T.  | PASS   | SVTYPE=BND | SGL                           |
| 3     | 5000 | m_s | G   | .G  | PASS   | SVTYPE=BND | SGL                           |
