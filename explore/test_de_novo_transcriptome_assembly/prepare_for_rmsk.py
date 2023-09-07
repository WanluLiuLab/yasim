
import json
import sys
from labw_utils.bioutils.parser.fasta import FastaIterator, FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.lwio.safe_io import get_writer


if __name__ == "__main__":
    basename = sys.argv[1]
    # with FastaIterator("sim/ce11_denovo_test_pbsim.fa") as fai:
    #     with FastaWriter("sim/ce11_denovo_test_pbsim_rmsk_prep.fa") as faw:
    #         for i, record in enumerate(fai):
    #             faw.write(
    #                 FastaRecord(
    #                     seq_id=str(i),
    #                     sequence=record.sequence
    #                 )
    #             )
    #             retl.append(record.seq_id)
    # with get_writer("sim/ce11_denovo_test_pbsim_rmsk_prep.json", is_binary=False) as w:
    #     json.dump(retl, w)
    retl = []

    with FastaIterator(f"{basename}.fa") as fai:
        with FastaWriter(f"{basename}_rmsk_prep.fa") as faw:
            for i, record in enumerate(fai):
                faw.write(
                    FastaRecord(
                        seq_id=str(i),
                        sequence=record.sequence
                    )
                )
                retl.append(record.seq_id)
    with get_writer(f"{basename}_rmsk_prep.json", is_binary=False) as w:
        json.dump(retl, w)
