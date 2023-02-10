# Copyright 2022 - Barcelona Supercomputing Center
# Author: Rodrigo Martin
# MIT License

from ..variants import VariantRecord

def _get_mate_id(variant_record: VariantRecord):
    if variant_record.alt_sv_breakend is None:
        raise ValueError('Variant record is not described in breakend notation')
    # Check if it has MATEID or PARID
    mate_id = variant_record.info.get('MATEID',
                                      variant_record.info.get('PARID',
                                                              f'{variant_record.alt_sv_breakend.contig}{variant_record.alt_sv_breakend.pos}'))
    mate_id = mate_id[0] if type(mate_id) != str else mate_id
    return mate_id


def _get_id(variant_record: VariantRecord):
    if 'MATEID' in variant_record.info or 'PARID' in variant_record.info:
        return variant_record.id
    else:
        return f'{variant_record.contig}{variant_record.pos}'


class PendingBreakends:
    def __init__(self):
        self.__pending_breakends = {}

    def push(self, variant_record: VariantRecord):
        alt_breakend_id = _get_mate_id(variant_record)
        breakend_id = _get_id(variant_record)
        # Check if alt is already in the dictionary
        previous_records = self.__pending_breakends.get(alt_breakend_id)
        if previous_records is None:
            # Create new entry for alt
            new_records = {breakend_id: variant_record}
            self.__pending_breakends[alt_breakend_id] = new_records
        else:
            # Already exists alt entry, add new entry to alt
            previous_records[breakend_id] = variant_record

    def pop(self, alt_variant_record: VariantRecord):
        breakend_id = _get_id(alt_variant_record)
        previous_records = self.__pending_breakends.get(breakend_id)
        if previous_records is None:
            return None
        previous_record_alt_breakend_id = _get_mate_id(alt_variant_record)
        previous_record_alt = previous_records.get(previous_record_alt_breakend_id)
        if previous_record_alt is None:
            return None
        previous_records.pop(previous_record_alt_breakend_id)
        if len(previous_records) == 0:
            self.__pending_breakends.pop(breakend_id)
        return previous_record_alt

    def remove(self, variant_record: VariantRecord):
        alt_breakend_id = _get_mate_id(variant_record)
        breakend_id = _get_id(variant_record)
        previous_records = self.__pending_breakends.get(alt_breakend_id)
        if previous_records is None:
            return
        previous_records.pop(breakend_id)
        if len(previous_records) == 0:
            self.__pending_breakends.pop(alt_breakend_id)

    def values(self):
        for alt_dict in self.__pending_breakends.values():
            for variant_record in alt_dict.values():
                yield variant_record

    def __len__(self):
        total_len = 0
        for alt_dicts in self.__pending_breakends.values():
            total_len += len(alt_dicts)
        return total_len
