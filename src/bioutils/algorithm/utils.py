import math
from typing import Tuple

coordinate_type = Tuple[float, float]


def euclid_distance(d1: coordinate_type, d2: coordinate_type) -> float:
    """
    >>> euclid_distance((0, 0), (3, 4))
    5.0
    """
    return math.sqrt((d1[0] - d2[0]) ** 2 + (d1[1] - d2[1]) ** 2)


def manhattan_distance(d1: coordinate_type, d2: coordinate_type) -> float:
    """
    >>> manhattan_distance((0, 0), (3, 4))
    7
    """
    return abs(d1[0] - d2[0]) + abs(d1[1] - d2[1])

