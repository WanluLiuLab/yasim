import os.path
from typing import Mapping, Any, Type

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor
from yasim.helper.depth import DepthType
from yasim.helper.llrg import AssembleDumb, pair_depth_info_with_transcriptome_fasta_filename, \
    generate_callback, AssembleSingleEnd, AssemblePairEnd
from yasim.llrg_adapter import BaseLLRGAdapter


def abstract_simulate(
        transcriptome_fasta_dir: str,
        output_fastq_prefix: str,
        depth: DepthType,
        jobs: int,
        simulator_name: str,
        adapter_args: Mapping[str, Any],
        assembler_args: Mapping[str, Any],
        adapter_class: Type[BaseLLRGAdapter],
        is_pair_end: bool,
        not_perform_assemble: bool = False
):
    output_fastq_dir = output_fastq_prefix + ".d"
    os.makedirs(output_fastq_dir, exist_ok=True)
    simulating_pool = ParallelJobExecutor(
        pool_name="Simulating jobs",
        pool_size=jobs
    )
    depth_info = list(pair_depth_info_with_transcriptome_fasta_filename(transcriptome_fasta_dir, depth))
    if not_perform_assemble:
        assembler = AssembleDumb(
            depth=depth,
            output_fastq_prefix=output_fastq_prefix,
            simulator_name=simulator_name,
            input_transcriptome_fasta_dir=transcriptome_fasta_dir,
            **assembler_args
        )
    elif is_pair_end:
        assembler = AssemblePairEnd(
            depth=depth,
            output_fastq_prefix=output_fastq_prefix,
            simulator_name=simulator_name,
            input_transcriptome_fasta_dir=transcriptome_fasta_dir,
            **assembler_args
        )
    else:
        assembler = AssembleSingleEnd(
            depth=depth,
            output_fastq_prefix=output_fastq_prefix,
            simulator_name=simulator_name,
            input_transcriptome_fasta_dir=transcriptome_fasta_dir,
            **assembler_args,
        )
    assembler.start()
    for transcript_depth, transcript_id, transcript_filename in tqdm(iterable=depth_info, desc="Submitting jobs..."):
        if transcript_depth == 0 or not os.path.exists(transcript_filename):
            continue
        sim_thread = adapter_class(
            input_fasta=transcript_filename,
            output_fastq_prefix=os.path.join(output_fastq_dir, transcript_id),
            depth=transcript_depth,
            **adapter_args
        )
        simulating_pool.append(sim_thread, callback=generate_callback(assembler, transcript_id))
    simulating_pool.start()
    simulating_pool.join()
    assembler.terminate()
    assembler.join()
