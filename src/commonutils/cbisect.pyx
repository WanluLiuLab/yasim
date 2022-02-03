# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: cbisect.pyx -- Cython implementation of bisect.
#
#  VERSION HISTORY:
#  2021-10-22 0.1  : Converted from standard Python library by YU Zhejian.
#  2021-10-23 0.2  : Clarified symantics and implemented by YUAN Ruihong.
#
# ==============================================================================

"""
cbisect.pyx -- Cython implementation of bisect.
"""
import bisect

cdef int cbisect_floor(int * a, int x, int lo, int hi):
    """
    Return ``pos`` such that ``a[pos] <= x < a[pos + 1]``
    
    :param a: A sorted array
    :param x: Value need to find
    :param lo: Lower bound
    :param hi: Higher bound
    :return: Left index of the value
    """
    if x < a[lo]:
        return lo - 1
    if x >= a[hi - 1]:
        return hi - 1

    cdef int mid
    while lo < hi - 1:
        # Invariant: a[lo] <= x < a[hi - 1]
        mid = (lo + hi) // 2
        # Invariant: lo < mid <= hi - 1
        # TODO: Benchmark trichotomy
        if a[mid] <= x:  # a[lo] < a[mid] <= x < a[hi - 1]
            lo = mid
        else:  # a[lo] <= x < a[mid] <= a[hi - 1]
            hi = mid
    return lo


cdef int cbisect_ceil(int * a, int x, int lo, int hi):
    """
    Return ``pos`` such that ``a[pos] < x <= a[pos + 1]``

    :param a: A sorted array
    :param x: Value need to find
    :param lo: Lower bound
    :param hi: Higher bound
    :return: Left index of the value
    """
    if x <= a[lo]:
        return lo
    if x > a[hi - 1]:
        return hi

    cdef int mid
    while lo < hi - 1:
        # Invariant: a[lo] < x <= a[hi - 1]
        mid = (lo + hi) // 2
        # Invariant: lo < mid <= hi - 1
        # TODO: Benchmark trichotomy
        if a[mid] < x:  # a[lo] < a[mid] < x <= a[hi - 1]
            lo = mid
        else:  # a[lo] < x <= a[mid] <= a[hi - 1]
            hi = mid
    return hi


def test_cbisect_left():
    cdef int in_arr[4]
    in_arr[:] = [1, 3, 5, 7]
    assert cbisect_floor(in_arr, 0, 0, 4) == -1
    assert cbisect_floor(in_arr, 1, 0, 4) == 0
    assert cbisect_floor(in_arr, 2, 0, 4) == 0
    assert cbisect_floor(in_arr, 3, 0, 4) == 1
    assert cbisect_floor(in_arr, 4, 0, 4) == 1
    assert cbisect_floor(in_arr, 5, 0, 4) == 2
    assert cbisect_floor(in_arr, 6, 0, 4) == 2
    assert cbisect_floor(in_arr, 7, 0, 4) == 3
    assert cbisect_floor(in_arr, 8, 0, 4) == 3

def test_cbisect_right():
    cdef int in_arr[4]
    in_arr[:] = [1, 3, 5, 7]
    assert cbisect_ceil(in_arr, 0, 0, 4) == 0
    assert cbisect_ceil(in_arr, 1, 0, 4) == 0
    assert cbisect_ceil(in_arr, 2, 0, 4) == 1
    assert cbisect_ceil(in_arr, 3, 0, 4) == 1
    assert cbisect_ceil(in_arr, 4, 0, 4) == 2
    assert cbisect_ceil(in_arr, 5, 0, 4) == 2
    assert cbisect_ceil(in_arr, 6, 0, 4) == 3
    assert cbisect_ceil(in_arr, 7, 0, 4) == 3
    assert cbisect_ceil(in_arr, 8, 0, 4) == 4
