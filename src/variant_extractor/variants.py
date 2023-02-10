# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# MIT License
from typing import NamedTuple, Optional, List, Dict, Any
from enum import Enum, auto


class VariantType(Enum):
    """Enumeration with the different types of variations
    """
    SNV = auto()
    DEL = auto()
    INS = auto()
    DUP = auto()
    INV = auto()
    CNV = auto()
    TRA = auto()
    SGL = auto()


class BreakendSVRecord(NamedTuple):
    """NamedTuple with the information of a breakend notated SV record
    """
    prefix: Optional[str]
    """Prefix of the SV record with breakend notation. For example, for :code:`G]17:198982]` the prefix will be :code:`G`"""
    bracket: str
    """Bracket of the SV record with breakend notation. For example, for :code:`G]17:198982]` the bracket will be :code:`]`"""
    contig: str
    """Contig of the SV record with breakend notation. For example, for :code:`G]17:198982]` the contig will be :code:`17`"""
    pos: int
    """Position of the SV record with breakend notation. For example, for :code:`G]17:198982]` the position will be :code:`198982`"""
    suffix: Optional[str]
    """Suffix of the SV record with breakend notation. For example, for :code:`G]17:198982]` the suffix will be :code:`None`"""


class ShorthandSVRecord(NamedTuple):
    """NamedTuple with the information of a shorthand SV record
    """
    type: str
    """One of the following, :code:`'DEL'`, :code:`'INS'`, :code:`'DUP'`, :code:`'INV'` or :code:`'CNV'`"""
    extra: Optional[List[str]]
    """Extra information of the SV. For example, for :code:`<DUP:TANDEM:AA>` the extra will be :code:`['TANDEM', 'AA']`"""


def _str_value(value):
    if isinstance(value, str):
        return value
    elif isinstance(value, float):
        return f'{value:.2f}'
    elif hasattr(value, '__iter__'):
        return ','.join([_str_value(v) for v in value])
    elif value is None:
        return '.'
    else:
        return str(value)


def _convert_info_key_value(key, value):
    if value is None:
        return key
    elif isinstance(value, bool):
        return key if value else None
    else:
        return key+'='+_str_value(value)


def _convert_sample_value(key, value):
    if key == 'GT':
        return '/'.join([_str_value(v) for v in value])
    else:
        return _str_value(value)


class VariantRecord(NamedTuple):
    """NamedTuple with the information of a variant record
    """
    contig: str
    """Contig name"""
    pos: int
    """Position of the variant in the contig"""
    end: int
    """End position of the variant in the contig (same as `pos` for TRA and SNV)"""
    length: int
    """Length of the variant"""
    id: Optional[str]
    """Record identifier"""
    ref: str
    """Reference sequence"""
    alt: str
    """Alternative sequence"""
    qual: Optional[float]
    """Quality score for the assertion made in ALT"""
    filter: List[str]
    """Filter status. PASS if this position has passed all filters. Otherwise, it contains the filters that failed"""
    info: Dict[str, Any]
    """Additional information"""
    format: List[str]
    """Specifies data types and order of the genotype information"""
    samples: Dict[str, Dict[str, Any]]
    """Genotype information for each sample"""
    variant_type: VariantType
    """Variant type"""
    alt_sv_breakend: Optional[BreakendSVRecord]
    """Breakend SV info, present only for SVs with breakend notation. For example, :code:`G]17:198982]`"""
    alt_sv_shorthand: Optional[ShorthandSVRecord]
    """Shorthand SV info, present only for SVs with shorthand notation. For example, :code:`<DUP:TANDEM>`"""

    def __str__(self):
        contig = self.contig
        pos = self.pos
        id_ = self.id if self.id else '.'
        ref = self.ref
        alt = self.alt
        qual = _str_value(self.qual)
        filter_ = ";".join(self.filter) if self.filter else '.'
        info_list = []
        for key, value in self.info.items():
            info_str = _convert_info_key_value(key, value)
            if info_str is None:
                continue
            info_list.append(info_str)
        if self.alt_sv_shorthand:
            info = info_list.insert(0, 'END='+str(self.end))
        info = ";".join(info_list)
        format_ = ":".join(self.format)
        samples_list = [":".join([_convert_sample_value(k, self.samples[sample_name][k])
                                 for k in self.format]) for sample_name in self.samples]
        samples = "\t".join(samples_list)
        return f'{contig}\t{pos}\t{id_}\t{ref}\t{alt}\t{qual}\t{filter_}\t{info}\t{format_}\t{samples}'.strip()
