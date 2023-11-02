"""
Transposon database that manageconserved transposon sequences.

.. versionadded:: 3.2.0
"""

import json
import random
import shutil

import h5py

from labw_utils.bioutils.parser.fasta import FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.lwio import get_writer
from labw_utils.commonutils.stdlib_helper import pickle_helper
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List, Mapping, Tuple, Optional

_lh = get_logger()


class TransposonDatabase:
    _accession_sequence_map: Mapping[str, str]
    _accession_hmm_map: Mapping[str, str]
    _accessions: List[str]
    _hmm_epool: List[Tuple[str, str]]

    @staticmethod
    def convert_dfam_hdf5(
        *,
        src_dfam_hdf5_file_path: str,
        required_txid: int,
        dst_index_file_path: str,
        dst_consensus_fa_path: Optional[str],
        dst_hmm_path: Optional[str],
        with_tqdm: bool = True,
        fetch_parent: bool = True,
    ) -> None:
        _lh.info("Enumerating taxons...")
        txids = [required_txid]
        accession_info = {}
        with h5py.File(src_dfam_hdf5_file_path) as ds:
            while fetch_parent:
                this_taxon_node = ds["Taxonomy"]["Nodes"][str(txids[-1])]
                if "Parent" in this_taxon_node:
                    txids.append(this_taxon_node["Parent"][0])
                else:
                    break

            _lh.info("Enumerating taxons finished with %d taxons in list", len(txids))
            _lh.info("Enumerating accessions...")
            accessions = set()
            for txid in txids:
                this_taxon_node = ds["Taxonomy"]["Nodes"][str(txid)]
                if "Families" in this_taxon_node:
                    accessions.update(this_taxon_node["Families"])
            if with_tqdm:
                it = tqdm(accessions)
            else:
                it = accessions
            for accession in it:
                try:
                    d = ds["Families"][accession[0:2]][accession[2:4]][accession[4:6]][
                        accession
                    ]
                except KeyError:
                    try:
                        d = ds["Families"]["ByName"][accession]
                    except KeyError:
                        _lh.warning(
                            "Accession '%s' neither resolved as name nor accession, skipped",
                            accession,
                        )
                        continue
                try:
                    accession_info[d.attrs["name"]] = {
                        "CONSENSUS": d.attrs["consensus"],
                        "HMM": d.attrs["model"],
                        "NAME": d.attrs.get("name", f"NAME-{accession}"),
                        "SUBTYPE": d.attrs.get(
                            "repeat_subtype", f"SUBTYPE-{accession}"
                        ),
                        "TYPE": d.attrs.get("repeat_type", f"TYPE-{accession}"),
                    }
                except KeyError:
                    _lh.warning(
                        "Name '%s' do not have consensus/model, skipped",
                        accession,
                    )
                    continue
            _lh.info(
                "Finished with %d accesions, writing...",
                len(accession_info),
            )
            pickle_helper.dump(
                accession_info,
                dst_index_file_path,
            )
            if dst_consensus_fa_path is not None:
                with FastaWriter(dst_consensus_fa_path) as faw:
                    for k, v in accession_info.items():
                        faw.write(FastaRecord(k, v["CONSENSUS"]))
            if dst_hmm_path is not None:
                with get_writer(dst_hmm_path, is_binary=False) as hmmw:
                    for v in accession_info.values():
                        hmmw.write(v["HMM"])
            _lh.info("Finished")

    @classmethod
    def load(cls, index_path: str, with_tqdm: bool = True):
        accession_info = pickle_helper.load(index_path, with_tqdm=with_tqdm)

        return cls(
            {k: v["CONSENSUS"] for k, v in accession_info.items()},
            {k: v["HMM"] for k, v in accession_info.items()},
        )

    def dump_json(self, dst_json_file_path: str):
        with get_writer(dst_json_file_path, is_binary=False) as w:
            json.dump(self._accession_sequence_map, w, indent=4)

    def __init__(
        self,
        accession_sequence_map: Mapping[str, str],
        accession_hmm_map: Mapping[str, str],
        **kwargs,
    ) -> None:
        _ = kwargs
        del kwargs
        self._accession_sequence_map = accession_sequence_map
        self._accession_hmm_map = accession_hmm_map
        self._accessions = list(self._accession_sequence_map.keys())
        self._hmm_epool = []

    def draw(self) -> Tuple[str, str]:
        """
        Randomly pick one transposon from the database.

        :returns: Transposon accession, start position and its sequence.
        """
        rdg = random.SystemRandom()
        selected_accn = rdg.choice(self._accessions)
        selected_seq = self._accession_sequence_map[selected_accn]
        return selected_accn, selected_seq
