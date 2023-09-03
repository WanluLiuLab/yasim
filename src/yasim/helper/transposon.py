"""
Transposon database that manageconserved transposon sequences.

.. versionadded:: 3.2.0
"""

import json
import random
import h5py 
from labw_utils.commonutils.lwio import get_writer
from labw_utils.commonutils.importer.tqdm_importer import tqdm

from labw_utils.typing_importer import Sequence, Tuple, Dict
from labw_utils.commonutils.stdlib_helper import pickle_helper 

class TransposonDatabase:
    _accession_sequence_map: Dict[str, str]

    @staticmethod
    def convert_dfam_hdf5(
        src_dfam_hdf5_path: str,
        required_txid: int,
        dst_index_path:  str, 
        with_tqdm:bool = True,
        fetch_parent: bool = True,
        ) -> None:
        txids = [required_txid]
        accession_sequence_map = {}
        ds = h5py.File(src_dfam_hdf5_path)
        while fetch_parent:
            this_taxon_node = ds["Taxonomy"]["Nodes"][str(txids[-1])]
            if "Parent" in this_taxon_node:
                txids.append(this_taxon_node["Parent"])
            else:
                break
        accessions = set()
        for txid in txids:
            accessions.update(ds["Taxonomy"]["Nodes"][str(txid)]["Families"])
        if with_tqdm:
            it = tqdm(accessions)
        else:
            it = accessions
        for accession in it:
            d = ds["Families"][accession[0:2]][accession[2:4]][accession[4:6]][accession]
            accession_sequence_map[d.attrs["name"]] = d.attrs["consensus"]

        pickle_helper.dump(accession_sequence_map, dst_index_path)

    @classmethod
    def load(cls, index_path: str, with_tqdm: bool = True):
        return cls(pickle_helper.load(index_path, with_tqdm=with_tqdm))

    def dump_json(self, dst_json_file_path: str):
        with get_writer(dst_json_file_path, is_binary = False) as w:
            json.dump(self._accession_sequence_map, w, indent = 4)

    def __init__(self, accession_sequence_map: Dict[str, str]) -> None:
        self._accession_sequence_map = accession_sequence_map

    def pick(self, min_length: int = 10) -> Tuple[str, int, str]:
        """
        Randomly pick one transposon from the database.

        :returns: Transposon accession, start position and its sequence.
        """
        rdg = random.SystemRandom()
        while True:
            selected_accn = rdg.choice(self._accession_sequence_map.keys())
            selected_seq = self._accession_sequence_map[selected_accn]
            selected_start = rdg.randint(0, len(selected_seq) - 1)
            selected_end = rdg.randint(selected_start, len(selected_seq))
            if selected_end - selected_start + 1 > min_length:
                return selected_accn, selected_start, selected_seq[selected_start:selected_end]
