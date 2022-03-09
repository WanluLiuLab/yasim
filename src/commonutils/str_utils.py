"""
str_utils.py -- String utilities

This file defines ANSI color supported by most terminals.
"""
import copy
import os

import string
from enum import Enum

from typing import Iterator, List, Dict, Any, Optional


class AnsiColorEnum(Enum):
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    PURPLE = "\033[35m"
    CRAYON = "\033[36m"
    CLEAR = "\033[0m"


class DumbColorEnum(Enum):
    RED = ""
    GREEN = ""
    YELLOW = ""
    BLUE = ""
    PURPLE = ""
    CRAYON = ""
    CLEAR = ""


def get_color(fd: int):
    """
    Get ANSI color dictionary for current file descriptor.

    Python can test whether the output is a tty. Other method have to employ ncurses.

    :param fd: File descriptor.
    :return: A color enum with a format like ``{RED = "\\033[31m"}``.
    """
    if os.isatty(fd):
        retd = AnsiColorEnum
    else:
        retd = DumbColorEnum
    return retd


def dict_exchange_key_val(in_dict: Dict[Any, Any]) -> Dict[Any, Any]:
    """
    To exchange the keys and values of one dictionary, that is,
    make a key value and make a value key.

    If the value contains duplicates, will use last key.

    >>> in_dict = {1: 'a', 2: 'b'}
    >>> dict_exchange_key_val(in_dict)
    {'a': 1, 'b': 2}

    >>> in_dict = {1: 'a', 2: 'a'}
    >>> dict_exchange_key_val(in_dict)
    {'a': 2}

    :param in_dict: Dictionary to be exchanged.
    :return: Exchanged dictionary.
    """
    out_dict = {}
    for k, v in in_dict.items():
        out_dict[v] = k
    return out_dict


def dict_translate(in_dict: Dict[str, Any], trans_dict: Dict[str, str]) -> Dict[str, Any]:
    """
    Dictionary Translator.

    This function will change the key of in_dict with the rules specified
    in trans_dict.

    For example:

    >>> dict_translate({'A':1, 'B':2, 'C':3}, {'A':'a', 'B':'b'})
    {'a': 1, 'b': 2, 'C':3}

    .. warning::
     The order of item will change!

    :param in_dict: The input dictionary.
    :param trans_dict: The translator.
    """
    new_dict = {}
    for old_key in in_dict.keys():
        if old_key in trans_dict.keys():
            new_key = trans_dict[old_key]
        else:
            new_key = old_key
        new_dict[new_key] = in_dict[old_key]
    return new_dict


def list_translate(in_list: List[str], trans_dict: Dict[str, str]) -> Iterator[str]:
    """
    List Translator.

    Translate the list as is specified in py:func:`dict_translate`.

    The order of the item will NOT be changed.

    >>> list_translate(['A', 'B', 'C'], {'A':'a', 'B':'b'})
    ['a', 'b', 'C']

    :param in_list: Input list
    :param trans_dict: The translator.
    :type trans_dict: dict
    :return: Translated dictionary.
    """
    for old_item in copy.deepcopy(in_list):
        if old_item in trans_dict.keys():
            yield trans_dict[old_item]
        else:
            yield old_item


def to_dict(
        in_str: str,
        field_sep: str = '\t',
        record_sep: str = '\n',
        quotation_mark: Optional[str] = None,
        resolve_str: bool = True
) -> Dict[str, Any]:
    """
    A simple parser to get key-value pairs to a dictionary.

    Key: string. Value: String, float or int.

    >>> in_str = '\\nCPU:\\t2\\nMEM:\\t5.1\\nPCIE:\\t3rd Gen\\nGRAPHICS:\\t"UHD630\\tRTX2070"\\nUSB: "3.1"\\nOthers:::info'
    >>> to_dict(in_str, field_sep=':', record_sep='\\n', quotation_mark="\\'\\"", resolve_str=True)
    {'CPU': 2, 'MEM': 5.1, 'PCIE': '3rd Gen', 'GRAPHICS': 'UHD630\\tRTX2070', 'USB': 3.1, 'Others': 'info'}
    >>> to_dict(in_str, field_sep=':', record_sep='\\n', quotation_mark="\\'\\"", resolve_str=False)
    {'CPU': '2', 'MEM': '5.1', 'PCIE': '3rd Gen', 'GRAPHICS': 'UHD630\\tRTX2070', 'USB': '3.1', 'Others': 'info'}
    >>> to_dict(in_str, field_sep=':', record_sep='\\n', quotation_mark=None, resolve_str=False)
    {'CPU': '2', 'MEM': '5.1', 'PCIE': '3rd Gen', 'GRAPHICS': '"UHD630\\tRTX2070"', 'USB': '"3.1"', 'Others': 'info'}
    >>> to_dict(in_str, field_sep=':', record_sep='\\n', quotation_mark=None, resolve_str=True)
    {'CPU': 2, 'MEM': 5.1, 'PCIE': '3rd Gen', 'GRAPHICS': '"UHD630\\tRTX2070"', 'USB': '"3.1"', 'Others': 'info'}

    :param in_str: Input string
    :param field_sep: Field separator, the FS variable in AWK programming language.
    :param record_sep: Record separator, the RS variable in AWK programming language.
    :param quotation_mark: If the key and value is not quoted, pass ``None``.
                           If quoted by single quote, pass ``'\''``.
                           If quoted by double quote, pass ``'\"'``.
                           If quoted by single and double quote, pass ``'\'\"'``.
                           Will not parse quoted by triple quotes.
    :param resolve_str: Whether to resolve strings to float or int.
    """
    retd = {}
    in_str_by_record = in_str.split(record_sep)
    for record in in_str_by_record:
        record = record.strip(string.whitespace + field_sep)
        lr = len(record)
        first_field_sep_pos = record.find(field_sep)
        if first_field_sep_pos == -1:
            continue
        record_key = record[0:first_field_sep_pos].rstrip()
        while not first_field_sep_pos == lr:
            if record[first_field_sep_pos] == field_sep:
                first_field_sep_pos += 1
            else:
                record_val = record[first_field_sep_pos:].lstrip()
                break
        else:
            record_val = ""
        if quotation_mark is not None:
            record_key = record_key.strip(quotation_mark)
            record_val = record_val.strip(quotation_mark)
        if resolve_str:
            try:
                if "." in record_val:
                    record_val = float(record_val)
                else:
                    record_val = int(record_val)
            except ValueError:
                pass
        retd[record_key] = record_val
    return retd
