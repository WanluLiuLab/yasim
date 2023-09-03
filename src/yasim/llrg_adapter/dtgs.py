"""
dtgs.py -- Wrapper of dTGS simulator (Dumb Third-Generation Sequencing Simulator).

.. versionadded:: 3.1.5
"""

__all__ = ("DTGSAdapter",)

from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.parser.fastq import FastqWriter
from labw_utils.bioutils.record.fastq import FastqRecord
from labw_utils.typing_importer import List, Union, Final, Mapping, Any
from yasim.llrg_adapter import BaseFunctionBasedLLRGAdapter


class DTGSAdapter(BaseFunctionBasedLLRGAdapter):
    """
    Wrapper of dTGS simulator (Dumb Third-Generation Sequencing Simulator).

    .. versionadded:: 3.1.5
    """

    def llrg_func(self):
        fav = FastaViewFactory(self.input_fasta, show_tqdm=False)
        with FastqWriter(self.output_fastq_prefix + ".fq") as fastq_writer:
            i = 0
            for chr_name in fav.chr_names:
                chr_seq = fav.sequence(chr_name)
                for _ in range(self.depth):
                    seqlen = len(chr_seq)
                    fastq_writer.write(FastqRecord(seq_id=str(i), sequence=chr_seq, quality="K" * seqlen))
                    i += 1

    llrg_name: Final[str] = "dtgs"
    _require_integer_depth: Final[bool] = True

    @staticmethod
    def validate_params(**kwargs) -> Mapping[str, Any]:
        return {}

    def __init__(
        self,
        *,
        src_fasta_file_path: str,
        dst_fastq_file_prefix: str,
        depth: Union[int, float],
        is_trusted: bool,
        other_args: List[str],
        **kwargs,
    ):
        """
        Initializer.

        :param src_fasta_file_path: Path of source FASTA.
        :param dst_fastq_file_prefix: Prefix od destination FASTQ.
        :param depth: Targeted sequencing depth. Would NOT be related to actual sequencing depth!
        :param other_args: Other arguments to be appended at the bottom of assembled CMD.
        :param is_trusted: Whether to skip input validation test.
        :raise LLRGInitializationException: On error.
        """
        _ = other_args, kwargs
        super().__init__(
            src_fasta_file_path=src_fasta_file_path,
            dst_fastq_file_prefix=dst_fastq_file_prefix,
            depth=depth,
            is_trusted=is_trusted,
            preserve_intermediate_files=None,  # Not effective.
        )

    def _pre_execution_hook(self) -> None:
        """Does not need extra preparation"""
        pass

    def _post_execution_hook(self):
        """Does not need extra preparation"""
        pass

    @property
    def is_pair_end(self) -> bool:
        return False
