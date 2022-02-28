import doctest
import functools
import math
from typing import Optional, List, Tuple

import numpy as np
import numpy.typing as npt


def smith_waterman_matrix(
        seq1: str, seq2: str,
        match_score: int = 5,
        mismatch_score: int = -4,
        indel_score: int = -4,
        is_global: bool = True
) -> npt.ArrayLike:
    l1 = len(seq1) + 1
    l2 = len(seq2) + 1
    score_matrix = np.zeros((l1, l2), dtype=int)
    for i in range(1, l1):
        for j in range(1, l2):
            if seq1[i - 1] == seq2[j - 1]:
                score = score_matrix[i - 1][j - 1] + match_score
            else:
                score = max(
                    score_matrix[i - 1][j] + indel_score,
                    score_matrix[i][j - 1] + indel_score,
                    score_matrix[i - 1][j - 1] + mismatch_score
                )
            if not is_global:
                score = max(score, 0)
            score_matrix[i][j] = score
    # warnings.warn(repr(score_matrix))
    return score_matrix


def smith_waterman_score(
        seq1: str, seq2: str,
        match_score: int = 5,
        mismatch_score: int = -4,
        indel_score: int = -4,
        is_global: bool = True
) -> int:
    """
    >>> smith_waterman_score('AAA', 'AAA')
    15
    >>> smith_waterman_score('AAA', 'ATA')
    6
    >>> smith_waterman_score('AAA', 'AA')
    10
    >>> smith_waterman_score('AA', 'TTTA')
    5
    """
    score_matrix = smith_waterman_matrix(seq1, seq2, match_score, mismatch_score, indel_score, is_global)

    return np.max(score_matrix[1:, 1:])


def smith_waterman_backtrack(
        seq1: str, seq2: str,
        match_score: int = 5,
        mismatch_score: int = -4,
        indel_score: int = -4,
        is_global: bool = True,
        seq1_title: str = "seq1",
        seq2_title: str = "seq2",
        alignment_title: str = "aln"
) -> Optional[List[str]]:
    """
    >>> smith_waterman_backtrack('AAA', 'AAA')[0]
    '>aln:seq1:qual:seq2:15\\nAAA\\nMMM\\nAAA'
    >>> smith_waterman_backtrack('AAA', 'ATA')[0]
    '>aln:seq1:qual:seq2:6\\nA--AA\\nMIIDD\\nATA--'
    >>> smith_waterman_backtrack('AAA', 'AA')[0]
    '>aln:seq1:qual:seq2:10\\nAAA\\nMMD\\nAA-'
    >>> smith_waterman_backtrack('AA', 'TTTA')[0]
    '>aln:seq1:qual:seq2:5\\n---AA\\nIIIMD\\nTTTA-'
    >>> smith_waterman_backtrack('JHBGTYBYTJAA', 'YJGYJTVJYAA')[0]
    '>aln:seq1:qual:seq2:10\\n-JHBGTYBY-T-J-AA\\nIMDDMDMDDIMIMIMM\\nYJ--G-Y--JTVJYAA'
    >>> smith_waterman_backtrack('TATATATGCGGGTAATTTAGGGCGGATCATGA', 'ATGCGGC')[0]
    '>aln:seq1:qual:seq2:30\\nTATATATGCGGGTAATTTAGGGCGGATCATGA\\nDMMDDDDMMMMMDDDDDDDDDDDDDDDDDDDD\\n-AT----GCGGC--------------------'
    """

    l1 = len(seq1)
    l2 = len(seq2)

    @functools.lru_cache()
    def location_inside_bound(input_location: Tuple[int, int]) -> bool:
        return 0 <= input_location[0] <= l1 and \
               0 <= input_location[1] and \
               input_location[1] <= l2

    score_matrix = smith_waterman_matrix(seq1, seq2, match_score, mismatch_score, indel_score, is_global)
    # warnings.warn(repr(score_matrix))
    best_match_locations = np.where(score_matrix == np.max(score_matrix[1:, 1:]))
    retl = set()
    for best_matchlocation in zip(*best_match_locations):
        if len(best_matchlocation) != 2:
            continue
        # Find a route to origin
        location = best_matchlocation
        step = []
        while location_inside_bound(location) and sum(location) != 0:
            step.append(location)
            next_score = -math.inf
            next_location = (-1, -1)
            for next_possible_location in [
                (location[0] - 1, location[1] - 1),
                (location[0], location[1] - 1),
                (location[0] - 1, location[1]),
            ]:
                if location_inside_bound(next_possible_location) and \
                        score_matrix[next_possible_location] > next_score:
                    next_score = score_matrix[next_possible_location]
                    next_location = next_possible_location
            location = next_location
        # Find a route to maximun point
        location = best_matchlocation
        step.reverse()
        while location_inside_bound(location) and sum(location) != l1 + l2:
            step.append(location)
            next_score = -math.inf
            next_location = (-1, -1)
            for next_possible_location in [
                (location[0] + 1, location[1] + 1),
                (location[0], location[1] + 1),
                (location[0] + 1, location[1]),
            ]:
                if location_inside_bound(next_possible_location) and \
                        score_matrix[next_possible_location] > next_score:
                    next_score = score_matrix[next_possible_location]
                    next_location = next_possible_location
            location = next_location
        step.append((l1, l2))
        # warnings.warn(repr(step))
        if len(step) <= 2:
            continue
        prev_location = (0, 0)
        out_array = []
        for this_location in step:
            if this_location == prev_location:
                continue
            elif this_location == (prev_location[0] + 1, prev_location[1] + 1):
                out_array.append((seq1[this_location[0] - 1], "M", seq2[this_location[1] - 1]))
            elif this_location == (prev_location[0] + 1, prev_location[1]):
                out_array.append((seq1[this_location[0] - 1], "D", "-"))
            elif this_location == (prev_location[0], prev_location[1] + 1):
                out_array.append(("-", "I", seq2[this_location[1] - 1]))
            else:
                out_array.append(("-", "?", "-"))
            prev_location = this_location
        rets = f">{alignment_title}:{seq1_title}:qual:{seq2_title}:{np.max(score_matrix[1:, 1:])}\n" + "\n".join(
            (seq1 + "" for seq1 in ["".join(seq) for seq in zip(*out_array)])
        )
        retl.add(rets)
    return list(retl)


def hamming_distance(str1: str, str2: str) -> int:
    """
    >>> hamming_distance("AAAA", "AATA")
    1
    """
    return sum(el1 != el2 for el1, el2 in zip(str1, str2))


def editing_distance(str1: str, str2: str) -> int:
    """
    >>> editing_distance("AAAA", "AAAA")
    0
    >>> editing_distance("AAAA", "AATA")
    1
    >>> editing_distance("AAAA", "AAA")
    1
    >>> editing_distance("AAAA", "TAAAA")
    1
    >>> editing_distance("A", "A")
    0
    >>> editing_distance("", "")
    0
    >>> editing_distance("", "A")
    1
    """
    l1 = len(str1)
    l2 = len(str2)
    if l1 == 0:
        return l2
    elif l2 == 0:
        return l1
    l1 += 1
    l2 += 1
    score_matrix = np.zeros((l1, l2), dtype=int)
    score_matrix[0, :] = range(l2)
    score_matrix[:, 0] = range(l1)
    for i in range(1, l1):
        for j in range(1, l2):
            if str1[i - 1] == str2[j - 1]:
                score_matrix[i][j] = score_matrix[i - 1][j - 1]
            else:
                score_matrix[i][j] = min(
                    score_matrix[i - 1][j - 1],
                    score_matrix[i][j - 1],
                    score_matrix[i - 1][j]
                ) + 1
    return score_matrix[l1 - 1][l2 - 1]


if __name__ == "__main__":
    doctest.testmod()
