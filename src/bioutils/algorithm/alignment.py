import functools
import math
from typing import Optional, List, Tuple

import numpy as np
import numpy.typing as npt


class SmithWatermanAligner:
    seq1: str
    """Sequence 1"""
    seq2: str
    """Sequence 2"""
    match_score: int
    """Score for match"""
    mismatch_score: int
    """Score for mismatch, should be negative"""
    indel_score: int
    """Score for insertions and deletions"""
    is_global: bool
    """Whether the alignment should be global"""

    _sw_matrix: npt.ArrayLike
    """The Smith-Waterman Alignment Matrix"""

    _score: Optional[int]
    """Final Smith-Waterman Alignment Score"""

    _backtrack: Optional[List[str]]
    """Final Smith-Waterman backtrack"""

    def __init__(self,
                 seq1: str, seq2: str,
                 match_score: int = 5,
                 mismatch_score: int = -4,
                 indel_score: int = -4,
                 is_global: bool = True
                 ):
        self.seq1 = seq1
        self.seq2 = seq2
        self.mismatch_score = mismatch_score
        self.match_score = match_score
        self.indel_score = indel_score
        self.is_global = is_global
        self._sw_matrix = None
        self._score = None
        self._backtrack = None

    def alignment(self) -> int:
        if self._sw_matrix is None:
            self._build_smith_waterman_matrix()
        self._score = np.max(self._sw_matrix[1:, 1:])
        return self._score

    @property
    def score(self) -> int:
        if self._score is None:
            return self.alignment()
        else:
            return self._score

    def _build_smith_waterman_matrix(self):
        l1 = len(self.seq1) + 1
        l2 = len(self.seq2) + 1
        self._sw_matrix = np.zeros((l1, l2), dtype=int)
        for i in range(1, l1):
            for j in range(1, l2):
                if self.seq1[i - 1] == self.seq2[j - 1]:
                    score = self._sw_matrix[i - 1][j - 1] + self.match_score
                else:
                    score = max(
                        self._sw_matrix[i - 1][j] + self.indel_score,
                        self._sw_matrix[i][j - 1] + self.indel_score,
                        self._sw_matrix[i - 1][j - 1] + self.mismatch_score
                    )
                if not self.is_global:
                    score = max(score, 0)
                self._sw_matrix[i][j] = score
        self._sw_matrix = self._sw_matrix

    def get_backtrack(self,
                      seq1_title: str = "seq1",
                      seq2_title: str = "seq2",
                      alignment_title: str = "aln"
                      ) -> Optional[List[str]]:
        """
        # >>> smith_waterman_backtrack('AAA', 'AAA')[0]
        # '>aln:seq1:qual:seq2:15\\nAAA\\n===\\nAAA'
        # >>> smith_waterman_backtrack('AAA', 'ATA')[0]
        # '>aln:seq1:qual:seq2:6\\nA--AA\\n=IIDD\\nATA--'
        # >>> smith_waterman_backtrack('AAA', 'AA')[0]
        # '>aln:seq1:qual:seq2:10\\nAAA\\n==D\\nAA-'
        # >>> smith_waterman_backtrack('AA', 'TTTA')[0]
        # '>aln:seq1:qual:seq2:5\\n---AA\\nIII=D\\nTTTA-'
        # >>> smith_waterman_backtrack('JHBGTYBYTJAA', 'YJGYJTVJYAA')[0]
        # '>aln:seq1:qual:seq2:10\\n-JHBGTYBY-T-J-AA\\nI=DD=D=DDI=I=I==\\nYJ--G-Y--JTVJYAA'
        # >>> smith_waterman_backtrack('TATATATGCGGGTAATTTAGGGCGGATCATGA', 'ATGCGGC')[0]
        # '>aln:seq1:qual:seq2:30\\nTATATATGCGGGTAATTTAGGGCGGATCATGA\\nD==DDDD====MDDDDDDDDDDDDDDDDDDDD\\n-AT----GCGGC--------------------'
        """
        if self._backtrack is not None:
            return self._backtrack
        if self._sw_matrix is None:
            self._build_smith_waterman_matrix()
        l1 = len(self.seq1)
        l2 = len(self.seq2)

        @functools.lru_cache()
        def location_inside_bound(input_location: Tuple[int, int]) -> bool:
            """Tell whether we're still inside the matrix"""
            return 0 <= input_location[0] <= l1 and 0 <= input_location[1] <= l2

        # warnings.warn(repr(self._sw_matrix))
        best_match_locations = np.where(self._sw_matrix == self.score)
        retl = set()
        for best_match_location in zip(*best_match_locations):
            # Assertion to ensure data correctness
            if len(best_match_location) != 2:
                continue
            # Find a route to origin
            location = best_match_location
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
                            self._sw_matrix[next_possible_location] > next_score:
                        next_score = self._sw_matrix[next_possible_location]
                        next_location = next_possible_location
                location = next_location
            # Find a route to maximum point
            location = best_match_location
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
                            self._sw_matrix[next_possible_location] > next_score:
                        next_score = self._sw_matrix[next_possible_location]
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
                    if self.seq1[this_location[0] - 1] == self.seq2[this_location[1] - 1]:
                        out_array.append((self.seq1[this_location[0] - 1], "=", self.seq2[this_location[1] - 1]))
                    else:
                        out_array.append((self.seq1[this_location[0] - 1], "M", self.seq2[this_location[1] - 1]))
                elif this_location == (prev_location[0] + 1, prev_location[1]):
                    out_array.append((self.seq1[this_location[0] - 1], "D", "-"))
                elif this_location == (prev_location[0], prev_location[1] + 1):
                    out_array.append(("-", "I", self.seq2[this_location[1] - 1]))
                else:
                    out_array.append(("-", "?", "-"))
                prev_location = this_location
            rets = f">{alignment_title}:{seq1_title}:qual:{seq2_title}:{np.max(self._sw_matrix[1:, 1:])}\n" + "\n".join(
                (_seq1 + "" for _seq1 in ["".join(_seq) for _seq in zip(*out_array)])
            )
            retl.add(rets)
        self._backtrack = list(retl)
        return self._backtrack


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
    seq1 = "TCTATTCTTTATGTCCCAAGAGGTTTTTATACGTTGTTGAAGCCACGATAGAAATAGCTCACTACGACTTCCTACTCGGTGGTTGCTCAGTTTTCAAGTCTTCTAAGCTCCTGGCCGAGTGTAACCGACCCTGGTCAGGTTCTGCTGTCGTGGCTGCAAGTCTATAGTAACGAAGTTGTATATCACTAAAGTATCTACTTCGTCTTATACGCTGTTTTGTCTACGACCCGGTAGAACCTTTGTGACTTAGGTAGCTTCTCTCATGTAACTTGATAAACTCTCGAGTAAGTTTCGGTATTTCCAGGCGTAACGGGTAAGCTACACTAAAAAATTTACTACTTTTTTTTGGTCTAACCTAACGCTAACCGCTTTTTGACTAATACCCGTCAGACTGTCCTTATGCCTGTAGACTTTTCGGTCCTTATTGGTGATTAGTCTACTGTCGCTAAGAAACTGCGAAAGGTAAAAGAAAAAACTACGGCAACGAACGTTAGGGCTCTTGAAATACTTTGCATAGTCTTATGAAATAGGTACTAACACTACTAACGCACTGGCCCCCTTGAAGATACTAGGGCTGCTGTCTCGTAGTTAACGTCATGTCGTGAATACCTTTCCTTATACGAACTAGTTCCTAAAGTACTCTCTCCGCAACAAGGGCCTTGAGCAAACTATGTTAAAGTTGTTGCCCACAGCCACGCAAGTCAGTTTGCTCGATTCCTACCTCACACAGGTCCGTTGAACGCAATAGAGCATACTACGCTAGGTGCTATGCTGGCTACCCAATAAGCTATTAAAGGACGGTGCAAAAATGTTTACGAGTTTCTAGGCTGCTACGTACAAGTCCTTCTCACGGGAAGTATAAAAAACCTTTTAGAACCTTTGTACAAAACCTGTCTTTAAGGGGCTTCTACAAGAACCAAGGGCTGGACAATTGTGAGAGCACGCTATACTGTACCCGCAAGCCTTCAACGTACGCTCTATAAACACTTTCACCGACTTACACGAAGTTTTAGTTCACTCTTAACTGTTAAAAACCATGTGGATTAAAGCACTCTAATAACTTACTAGGCTTGTTTACCCACAGCAAGTAACGAATAGTTGTTGAACCTGGTTAATGTAAAGTGTAAACGACTAAGGCCTATAATGCCCTGAACTTGAAGCTTACCTGTTCAGCAAGCTGTCTCTAGGTCATCATAATCTTGGTAGTCATTGACTTGTCCTTCGTAGGTATTTAAGATGTAAAAGCCTCCTGTAGCCTACTAAAATCTTGTGAAGGAGTTTCTATGCACCGACGGACTGTCAACCTAAAATTACTACTGAAACTTATAGTTAACCCTAGTTAGCTATTATGTGTCTTCAGTTTTCCTTAACTGTTGACACCTCCTGTCGTACACCACTTTTGAGGCTTGCTCCCCTTGCTCAATACCTTTCCCCAAGACAGACCTTAGTTTTAAAACCTAGCAACCTTGCTTCGAGCGGCTTTAGTTTTCTTTCAGCTTTTTTTAACGACGGTGGTGAAAAGTTCAGGTAGCTTTGCTAGTGGTGATAGTGAAGTTCACTCAGTCAGACTAACTATTTAAGTCGTTTTCTGTTAAGGCTTGTTAGGTGGTCGCTAAGGGTCATTAAACTAAGACGGTTGAGCCAAGGTGTAATAGTAGATTAAGGGGAGATTGTTTATTGTCGAGTCTCGTCGTCTATGGCTGCCCTTAGCGACGTATTCACGTCAGTTGCTACTAAACTGATGTTAATCCTGAACTTCTGCTACTGAGTGGTCTACGTCGATTCTAAGGGCTACTAAAGTAGTGCCTTTGTCCCTGCTGTTCAGTGCACCAGCCGCATAACGGTGGTTGCGCACCAGCTCGAACTATATAGCCGTCTTCTATAAAGAGGCCGCTGAACCTCCGGAGATGGTGGCTCCGAAAAAGGGCTTCTTTTTCCCGCTGTGTGTCGTCGTCGTCGTCAGGCTCTCTTAGTGGTTGTTTAGCGGTGACGACAGCAGCGCTTTTTTTCCCCTGAAAAAAACTAAAAAGCCGCGATAACTTTAGTTTAGTAGCAGAAGAACTAGTAAGCCTCTTTAAAGGTCGGCCTGTCGGTTAGTTTGTCGTTTTTTTCGTCTACTACTATTGCTTATGGCTTACTTGCATTGTTATGACCTTGTAGACCAGGTTTCTATTGTGACTACTTCTAGTTGCTTACTAAGTCTATCACTGCTACCTCCTATTACTAAAACTTTTCAGAAAACAACGACTGCTTCTACTGCCTCTTGTTCTCTATCGGTAAAGTAACTTTATTTAGTAGACGAAGAGTCTGAAACTTGAGTCAAAGCCGTCGGTGAGGTTGTAGAGTCCTCGGGTTCCATAAGATTACTCAAAGATAAACTGAATAACCAGCCAGTTGTTGATACTATAAGCAGGATAGCTTCTACTCGAGTTTGGAAAACCTAGGAGTGAGGAGTAAAGAGGTTAAGAGGTTAAGATTGTCAATAAACCTCTCGGCTAAGCTACACTACTTCCAATTAGATAGCGGTTTCTCGAGCTAGCTCACAGTAAACCTTGATTACCAGAGGCTATAGGGGAAGTGGTGCTTCTTTTCGGTTGTCATGTAGACTGGTAGAGGTGGTTTGTGGCCTTGAACTTTTCGCTATAGAGGGAACATTACGAGGTGTATGCCTTCTTTTTCCTTATTGACCCTTAAAATCTAAACTACTTTGCTTGTACAGTTTAATGCTTAAGTAGCTAAGTACCTCTTCGGCGTTACCTGCCCTTCTAACGGCCTATAATCAGATTGAAGCTAACTTTGTCCGTTCGTGAGGGTTGTACAAGGACAAATTAAAGAACAACTTAGAATATACTCTCAGTACCCCGGCTAACCTCCAGTAACCTTCGCCTCCCCACCATTACCTTCAAGATGGTTGGTATCAATGACCTCTAGGTAGCGAACGCCCCTCTGGATAGGGCAGTAACACGAGTGTTGAAAGGGCGTTCACACGATATGCGAACCGCACGACGTACCGG"
    seq2 = "AGTAGAGCAGTGCGAAAAGACCGCTGTCCTTGAGTTTAGTATCTTTTGACCGTTGTTTAGAAGTAGACACTAGAGGACAGGTTCGTATAAACGATAAAGTACTTAGCTTACTCTATAGCTGACTGCTTGAGCTCGGACGAGTAGGTAGCGTTCGAAGGCCTTATGTTCGAATTCGGCTCCGGTACGTCGTGCAGGTTCGCATATCGTGTGACGCCCTTTCAACACTCGTGTTACTGCTATCCAGAGGCGCGTTTCGCTACCTAGAGGTCATCGATACCAACATCTTGAAGGTAATGGTGGTGGAGGCGAAGGTGATGGAGGTTAGCCGGGGTACTGAGAGTATATTCTAAGTTGTTCTAGATGACCTTGTACAACCTCACGACGGACAAGTTAGCTTCAATCTGATTATAGGCTTAGAAGGCAGTAACGCCGAAGAGGTACTTAGCTACTTAAGCTACTATAACTGTACAAGTAAAGTAGCTTTAGATTTTAAGGGTCAAAGAGGAAAAGAAGGCATAACTCGTAATGTTACCTCAGTAGCGAAAAGTTCAAAGGCCACAAACCACCTCTACCAAGTCTACATGACAACCGAAAAGAAGCACCACTTCACCTATAGCCTCTAGGTAATCAAGGTTTACTGTGACTAGCTCGAGAACCGCTTTCTAATTGGAAGTAGTGTAGCTTAGCCGAGAGGTTTATTGACAATCTTCAACCTCTTAACCATTTACTCCTACTCCTAGGTTTTCCAAACTCGAGTAGAAGCTATCTGCTTTATCGTATCAAACTGCTGGTTATTCAGTTACTTTTGAGTAATCTATGACCGAGGACTCTACAACCTCACCGACGGCTTTGACTAAGTTTCAGACTGCTTCTTCTACTAAATTAAAGTTAACTATTACCTGATAGCAGAACAAGAGGCAGTAGAAGCAGTCGTTGTTTTCTGAAAAGTTTTAGTAATAGGAGTAGCAGTGATAGACTTAGTAAGCAACTTAGAAGTAGTCACAACTTGAAACCTGGTCTCAAAAGGTACATAACGACGGCAAGTAGCCATAAGTCGAATAGTAGTAGACGAAAAAGAACGACAACTAACCGACAGGCCTGACCTTAACAGAGGCTTACTAGTCTTCTGCTACTAAACTAAAGTTATCCGCGGCGTTTTAGTTTGTTTTCAGGGAAAAAAGCGCTGCTGTCGTCACCGCTAGAACAACCACTAAAAGAGCCTGACGACGACGACGACACACAGCGGAAAAAGAAGCCTTTCGAGCCACCATCTCGAGGTTCAGCGGCCTCTTTATAGAAGACGGCTTATAGTTCGAAGCTGGTGCGCAACCACCGTTATGCGGCTGGTGCACCTGAACAGCAGGCACAAAGGCTACACTTTAGTAGCCTTAGAAGTCGACGTAGACCACTCAGTAGCAGAAGTTCAGGATTAAATCAGTTTAGTAGTCAACTTGACGTGAAACGTCGCTAAGGGCAGCTTAAGACGACGAGACTCGACAATAACAATCTCCTTAATCTACTATTACACCTTAGGCTCAACCGTCTTAGTTTATGACCCTTAGCGACACCTAACAAGCCTAACGAAGAAAACGACTTAAATAGTTAGTCTGACTGACGTGACACTTCACTATCACCACTAGCAAAGCTACCTGAACTTTTACACCGTCGTTAAAAAAGCTGAAGAAAACTAAAGCCGCTCGAAGCAAGGTTGCTAGGTTTGTAAAACAACTAAGGTCTGTCTTGCGACAGGTATTGAGCAAGGGGAGCAAAGCCTCAAAAGTGGTTACGACAGCAGTGTCAACAGTTAAGGAAAACTAAGACACATAATAGCTAACTAGGGTTAACTATAAGTTTCAGTAGTAATTTTAGTTGAAGTCCGTCGGTGCATAGAAATCCTTCACAAGATTTTAGTAGGCTACACAGAAGGCTTTTACATCTTAAACTACCTAGAAGGACAAGTCAATGACTACCAAGATTATGATGACCTAGAGACTGTTGCTGAACAGGTAAAGCTTCAAGTTCAGGGCATTATAGGCCTAATCGTTTATCACTTTACTTACCAGGTTCAACAACTATTCGTTACTTGCTGTGGGTAAACAAGCCTAGTAAGTTATTAGAGTGCTTTAATCCACTGTTTTACAGTTAAGAGTGAACTAAACTTCGTGTAAGTCGGTGAAAGTGTTTATAGAGCGTACGTTGAAAGGCTTGGCGGGTACAGTATAGCGTGCTCTCACATTAGTCCACCCCTTGGTTCTTGTGAAGCCACTTAAAGACAGGTTTGTACAAAGGTTCTAAAGGTTTTTATACTTCCCGTGAGAAGGACTGTACGTAGCAGCCTAGAAACTCGTAAACACTTTTTGCACCGTCCTTTAATAGCTTATTGGTAGCAGCTAGCACCTAGCGTAGTATGCTCTATTGGACCGTTTCACGGACCTGTGTGAGGTAGGAATCTAGCAAACTGACCTTGCGTGGCTGTGGACTAACTATTTACAATTTCTCAGGCCCTTGTTGCGGTAGAGAGTACTTTAGGACTAGTTCGTATAAGGAAAGGTATTCACGGACATGACGTAACTACGAGACAGCAGCCTAGTATCTCAAAGGCCAGTAGCGTTAGTAGTGTTAGTACCTATTTCATAAGACTATAAGTATTTCAAGGCCTAACGTTCGTTGCCGTAGTTTTTCTTTACCTTATCGCAGTTTATAGCGACAGTAGACTAATCACCAATAAGGACCGAAAAGTCTACAACATAAGTGACAGTCTGACGGGTATTAGTCAAAAGCGGCTTAGCGTAGGGTTAGACCAAAAAAAAGTAGTAAATTTTAGTGTAGCTTACCGTTAGCCTGGAGATACCGAAACTTACTCGAGAGTTTATCAAGTTACATGAGAAGCTATCCTAAGTCACAAAGGTTCTACCGGTCGTAGACAAAACAGCGTATAAGACGAAGTAGATACTTTAGTGATATACAACTTCGTTACTATAGACCTTCGCAGCACGACAGAGAACCTGACAGGTCGGTTACACCTCGGCCAGAGCTTAGAAGACTTGAAAACTGAGAACCACCGAGTAGAAGTCGTAGTAGCTATTTCTATCGTGGCGGTTCAACAACGTTTAAACTCTTGGGACATAAAGAAAATAGGTA"

    a = SmithWatermanAligner(
        seq1, seq2,
        mismatch_score=-10,
        indel_score=-50,
        is_global=False
    )
    # print(a.get_backtrack()[0])

    b = SmithWatermanAligner(
        'JHBGTYBYTJAA', 'YJGYJTVJYAA',
        # mismatch_score=-10,
        # indel_score=-50,
        # is_global=False
    )
    print(b.get_backtrack()[0])
    import matplotlib.pyplot as plt
    plt.imshow(b._sw_matrix)
    plt.show()
