# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_str_utils.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-26 0.1  : Added by YU Zhejian.
#
# ==============================================================================
"""
test_str_utils.py -- Unit test of corresponding module.
"""
from commonutils import str_utils


def test_dict_exchange_key_value():
    in_dict = {"a": "1", 2: "b"}
    out_dict = str_utils.dict_exchange_key_val(in_dict)
    assert out_dict == {"b": 2, "1": "a"}


def test_dict_translate():
    in_dict = {"A": 1, "B": 2, "C": 3}
    trans_dict = {"A": "a", "B": "b"}
    out_dict = str_utils.dict_translate(in_dict, trans_dict)
    assert out_dict == {"a": 1, "b": 2, "C": 3}


def test_list_translate():
    in_list = ["A", "B", "C"]
    trans_dict = {"A": "a", "B": "b"}
    out_list = list(str_utils.list_translate(in_list, trans_dict))
    assert out_list == ["a", "b", "C"]


def test_to_dict():
    in_str = """
    CPU:	2
    MEM:	5.1
    PCIE:	3rd Gen
    GRAPHICS:	UHD630	RTX2070
    """
    out_dict = str_utils.to_dict(in_str, field_sep=':', record_sep='\n')
    assert out_dict == {"CPU": 2, "MEM": 5.1, "PCIE": "3rd Gen", "GRAPHICS": "UHD630\tRTX2070"}
