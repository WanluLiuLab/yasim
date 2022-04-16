import bisect

from bioutils.datastructure._gv_errors import *
from bioutils.datastructure.gene_view_proxy import DEFAULT_SORT_EXON_EXON_STRAND_POLICY, Transcript, Exon, Gene


class TranscriptMutator:

    @staticmethod
    def del_exon(transcript: Transcript, exon_number:int):
        transcript._exons.pop(exon_number)

    @staticmethod
    def sort_exons(
            transcript: Transcript,
            exon_number_policy: str = DEFAULT_SORT_EXON_EXON_STRAND_POLICY,
            remove_exon_duplicates: bool = True
    ):
        if transcript.number_of_exons == 0:
            return
        transcript._exons = sorted(transcript._exons)
        if remove_exon_duplicates and transcript.number_of_exons >= 2:
            exon_id = 1
            while exon_id < len(transcript._exons):
                if transcript.get_nth_exon(exon_id) == transcript.get_nth_exon(exon_id - 1):
                    transcript._exons.pop(exon_id)
                else:
                    exon_id += 1
        if exon_number_policy == "stranded":
            if transcript.strand == '+':
                for i in range(len(transcript._exons)):
                    transcript._exons[i].exon_number = i + 1
            elif transcript.strand == '-':
                for i in range(len(transcript._exons)):
                    transcript._exons[len(transcript._exons) - i - 1].exon_number = i + 1
        elif exon_number_policy == "unstranded":
            for i in range(len(transcript._exons)):
                transcript._exons[i].exon_number = i + 1

    @staticmethod
    def fast_add_exon(
            transcript: Transcript,
            exon: Exon
    ):
        TranscriptMutator.add_exon(transcript, exon, False, False, False)

    @staticmethod
    def add_exon(
            transcript: Transcript,
            exon: Exon,
            check_duplicate: bool = True,
            check_same_chrome: bool = True,
            check_same_strand: bool = True
    ):
        new_pos = bisect.bisect_left(transcript._exons, exon)
        if check_duplicate and new_pos < len(transcript._exons) and exon == transcript._exons[new_pos]:
            raise DuplicatedExonError
        if check_same_chrome and exon.seqname != transcript.seqname:
            raise ExonInATranscriptOnDifferentChromosomeError
        if check_same_strand and exon.strand != transcript.strand and exon.strand != "." and transcript.strand != ".":
            raise ExonInATranscriptOnDifferentStrandError
        transcript._exons.insert(new_pos, exon)


class GeneMutator:

    @staticmethod
    def del_transcript(gene: Gene, transcript_id: str):
        index = gene._transcript_ids.index(transcript_id)
        gene._transcript_ids.pop(index)
        gene._transcripts.pop(index)

    @staticmethod
    def replace_transcript(
            gene: Gene,
            transcript: Transcript,
            **kwargs
    ):
        gene.del_transcript(transcript.transcript_id)
        gene.add_transcript(transcript, **kwargs)

    @staticmethod
    def fast_add_transcript(
            gene: Gene,
            transcript: Transcript,
    ):
        GeneMutator.add_transcript(gene, transcript, False, False, False)

    @staticmethod
    def add_transcript(
            gene: Gene,
            transcript: Transcript,
            check_duplicate: bool = True,
            check_same_chrome: bool = True,
            check_same_strand: bool = True
    ):
        if transcript.transcript_id in gene._transcript_ids:
            raise DuplicatedTranscriptIDError
        new_pos = bisect.bisect_left(gene._transcripts, transcript)
        if check_duplicate and new_pos < len(gene._transcripts) and transcript == gene._transcripts[new_pos]:
            raise DuplicatedTranscriptError
        if check_same_chrome and transcript.seqname != gene.seqname:
            raise TranscriptInAGeneOnDifferentChromosomeError
        if check_same_strand and transcript.strand != gene.strand and transcript.strand != "." and gene.strand != ".":
            raise TranscriptInAGeneOnDifferentStrandError
        gene._transcript_ids.insert(new_pos, transcript.transcript_id)
        gene._transcripts.insert(new_pos, transcript)
