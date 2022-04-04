# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martín Posada
# BSC AS IS License
from typing import NamedTuple, Optional, List
from enum import Enum, auto
import pysam


class VariantType(Enum):
    """Enumeration with the different types of variations
    """
    SNV = auto()
    DEL = auto()
    INS = auto()
    DUP = auto()
    INV = auto()
    CNV = auto()
    TRN = auto()
    SGL = auto()


class BracketSVRecord(NamedTuple):
    """NamedTuple with the information of a bracketed SV record
    """
    prefix: Optional[str]
    """Prefix of the SV record with bracket notation. For example, for :code:`G]17:198982]` the prefix will be :code:`G`"""
    bracket: str
    """Bracket of the SV record with bracket notation. For example, for :code:`G]17:198982]` the bracket will be :code:`]`"""
    contig: str
    """Contig of the SV record with bracket notation. For example, for :code:`G]17:198982]` the contig will be :code:`17`"""
    pos: int
    """Position of the SV record with bracket notation. For example, for :code:`G]17:198982]` the position will be :code:`198982`"""
    suffix: Optional[str]
    """Suffix of the SV record with bracket notation. For example, for :code:`G]17:198982]` the suffix will be :code:`None`"""


class ShorthandSVRecord(NamedTuple):
    """NamedTuple with the information of a shorthand SV record
    """
    type: str
    """One of the following, :code:`'DEL'`, :code:`'INS'`, :code:`'DUP'`, :code:`'INV'` or :code:`'CNV'`"""
    extra: List[str]
    """Extra information of the SV. For example, for :code:`<DUP:TANDEM:AA>` the extra will be :code:`['TANDEM', 'AA']`"""


class VariantRecord(NamedTuple):
    """NamedTuple with the information of a variant record
    """
    contig: str
    """Contig name"""
    pos: int
    """Position of the variant in the contig"""
    end: int
    """End position of the variant in the contig (same as `pos` for TRN and SNV)"""
    length: int
    """Length of the variant"""
    id: str
    """Record identifier"""
    ref: str
    """Reference sequence"""
    alt: str
    """Alternative sequence"""
    qual: Optional[float]
    """Quality score for the assertion made in ALT"""
    filter: pysam.VariantRecordFilter
    """Record filter"""
    info: pysam.VariantRecordInfo
    """Dictionary of information fields"""
    variant_type: VariantType
    """Variant type"""
    alt_sv_bracket: Optional[BracketSVRecord]
    """Bracketed SV info, present only for SVs with bracket notation. For example, :code:`G]17:198982]`"""
    alt_sv_shorthand: Optional[ShorthandSVRecord]
    """Shorthand SV info, present only for SVs with shorthand notation. For example, :code:`<DUP:TANDEM>`"""
