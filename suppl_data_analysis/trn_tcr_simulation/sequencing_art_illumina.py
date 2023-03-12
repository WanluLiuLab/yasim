
from yasim.llrg_adapter.art import ArtAdapter

if __name__ == "__main__":
    sim_thread = ArtAdapter(
        input_fasta="sim.fa",
        output_fastq_prefix="sim",
        depth=40.0,
        exename="art_illumina",
        is_pair_end=False,
        mflen_std=0,
        mflen_mean=0,
        other_args=[],
        rlen=150,
        sequencer="HS25"
    )
    sim_thread.start()
    sim_thread.join()

