"""
Transposon database that manageconserved transposon sequences.

.. versionadded:: 3.2.0
"""

import json
import random
import h5py
from labw_utils.commonutils.lwio import get_writer
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

from labw_utils.typing_importer import List, Tuple, Dict
from labw_utils.commonutils.stdlib_helper import pickle_helper

_lh = get_logger()


class TransposonDatabase:
    _accession_sequence_map: Dict[str, str]
    _accessions: List[str]

    @staticmethod
    def convert_dfam_hdf5(
        src_dfam_hdf5_file_path: str,
        required_txid: int,
        dst_index_file_path: str,
        with_tqdm: bool = True,
        fetch_parent: bool = True,
    ) -> None:
        _lh.info("Enumerating taxons...")
        txids = [required_txid]
        accession_sequence_map = {}
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
                    d = ds["Families"][accession[0:2]][accession[2:4]][accession[4:6]][accession]
                except KeyError:
                    try:
                        d = ds["Families"]["ByName"][accession]
                    except KeyError:
                        _lh.warning("Accession '%s' neither resolved as name nor accession, skipped", accession)
                        continue
                accession_sequence_map[d.attrs["name"]] = d.attrs["consensus"]

            _lh.info("Finished with %d accesions, writing...", len(accession_sequence_map))

            pickle_helper.dump(accession_sequence_map, dst_index_file_path)
            _lh.info("Finished")

    @classmethod
    def load(cls, index_path: str, with_tqdm: bool = True):
        return cls(pickle_helper.load(index_path, with_tqdm=with_tqdm))

    def dump_json(self, dst_json_file_path: str):
        with get_writer(dst_json_file_path, is_binary=False) as w:
            json.dump(self._accession_sequence_map, w, indent=4)

    def __init__(self, accession_sequence_map: Dict[str, str]) -> None:
        self._accession_sequence_map = accession_sequence_map
        self._accessions = list(self._accession_sequence_map.keys())

    def draw(self) -> Tuple[str, str]:
        """
        Randomly pick one transposon from the database.

        :returns: Transposon accession, start position and its sequence.
        """
        rdg = random.SystemRandom()
        selected_accn = rdg.choice(self._accessions)
        selected_seq = self._accession_sequence_map[selected_accn]
        return selected_accn, selected_seq
