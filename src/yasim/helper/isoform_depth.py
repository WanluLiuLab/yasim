"""
isoform_depth.py -- Miscellaneous isoform depth functions.

.. versionadded:: 3.2.1
"""

__all__ = ("generate_isoform_depth",)

import numpy as np

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe
from labw_utils.typing_importer import Sequence, Mapping
from yasim.helper import depth

_lh = get_logger(__name__)


def generate_isoform_depth(
    gene_id_to_transcript_ids_map: Mapping[str, Sequence[str]],
    gene_level_depth: depth.DepthType,
    alpha: int = depth.DEFAULT_ALPHA,
    low_cutoff: float = depth.DEFAULT_LOW_CUTOFF,
    high_cutoff_ratio: float = depth.DEFAULT_MU,
) -> depth.DepthType:

    transcript_level_depth = {}
    for gene_id, transcript_ids in gene_id_to_transcript_ids_map.items():
        if gene_id not in gene_level_depth:
            _lh.warning(
                "GEN ISOFORM DEPTH: Gene %s defined in GTF but not gene-level depth",
                gene_id,
            )
        if gene_level_depth[gene_id] == 0:
            for transcript_id in transcript_ids:
                transcript_level_depth[transcript_id] = 0
            continue
        try:
            this_transcript_level_depth = depth.simulate_isoform_variance_inside_a_gene(
                n=len(transcript_ids),
                mu=gene_level_depth[gene_id],
                low_cutoff=low_cutoff,
                alpha=alpha,
                high_cutoff_ratio=high_cutoff_ratio,
            )
        except depth.GenerationFailureException:
            _lh.error(
                "GEN ISOFORM DEPTH: Generation failed for gene %s -- SKIPPED",
                gene_id,
            )
            continue
        for i, transcript_id in enumerate(transcript_ids):
            transcript_level_depth[transcript_id] = this_transcript_level_depth[i]
    _lh.info(
        "GEN ISOFORM DEPTH: Generation of isoform-level depth: Final distribution: %s",
        describe(np.array(list(transcript_level_depth.values()))),
    )
    return transcript_level_depth
