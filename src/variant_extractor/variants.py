# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# MIT License
from functools import cached_property
from typing import NamedTuple, Optional, List, Dict, Any
from enum import Enum, auto
import copy

import pysam


def _build_filter(rec: pysam.VariantRecord) -> List[str | int]:
    return [f for f in rec.filter]


def _build_info(rec: pysam.VariantRecord) -> Dict[str, Any]:
    info = dict()
    for key, value in rec.info.items():
        info[key] = value
    return info


def _build_format(rec: pysam.VariantRecord) -> List[str]:
    return [f for f in rec.format]


def _build_samples(rec: pysam.VariantRecord) -> Dict[str, Dict[str, Any]]:
    samples = dict()
    for sample_name in rec.samples:
        sample_dict = dict()
        for key, value in rec.samples[sample_name].items():
            sample_dict[key] = value
        samples[sample_name] = sample_dict
    return samples


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


class VariantRecord():
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
    filter: List[str | int]
    """Filter status. PASS if this position has passed all filters. Otherwise, it contains the filters that failed"""
    variant_type: VariantType
    """Variant type"""
    alt_sv_breakend: Optional[BreakendSVRecord]
    """Breakend SV info, present only for SVs with breakend notation. For example, :code:`G]17:198982]`"""
    alt_sv_shorthand: Optional[ShorthandSVRecord]
    """Shorthand SV info, present only for SVs with shorthand notation. For example, :code:`<DUP:TANDEM>`"""

    def __init__(self, rec: pysam.VariantRecord, contig: str, pos: int, end: int,
                 length: int, id: Optional[str], ref: str,
                 alt: str, variant_type: VariantType,
                 alt_sv_breakend: Optional[BreakendSVRecord] = None,
                 alt_sv_shorthand: Optional[ShorthandSVRecord] = None):
        self._rec = rec
        self.contig = contig
        self.pos = pos
        self.end = end
        self.length = length
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = rec.qual
        self.filter = _build_filter(rec)
        self.variant_type = variant_type
        self.alt_sv_breakend = alt_sv_breakend
        self.alt_sv_shorthand = alt_sv_shorthand

    @cached_property
    def info(self):
        """Additional information"""
        return _build_info(self._rec)

    @cached_property
    def format(self):
        """Specifies data types and order of the genotype information"""
        return _build_format(self._rec)

    @cached_property
    def samples(self):
        """Genotype information for each sample"""
        return _build_samples(self._rec)

    def _replace(self, **kwargs):
        self_copy = copy.copy(self)
        for key, value in kwargs.items():
            setattr(self_copy, key, value)
        return self_copy

    def __str__(self):
        contig = self.contig
        pos = self.pos
        id_ = self.id if self.id else '.'
        ref = self.ref
        alt = self.alt
        qual = _str_value(self.qual)
        filter_ = ";".join(map(str, self.filter)) if self.filter else '.'
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
