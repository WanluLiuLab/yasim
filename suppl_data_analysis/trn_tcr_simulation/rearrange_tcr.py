"""

::
    [    V    ][    J    ]
    [V     ]      [     J]
        XXXX      XXX
        [    CDR3   ]
        C           F
           [ INS ]


[INS]: Cdr3InsertionTable
[X]: Cdr3DeletionTable
"""

import itertools
import json
import os
import random
import uuid
from typing import Dict, List, Tuple

from create_tcr_cache import TCRTranslationTableType
from labw_utils.bioutils.algorithm.sequence import TRANSL_TABLES, TRANSL_TABLES_NT
from labw_utils.commonutils.appender import load_table_appender_class
from labw_utils.commonutils.appender.typing import TableAppenderConfig
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_reader, get_writer
from labw_utils.commonutils.shell_utils import rm_rf
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

Cdr3DeletionTableType = Dict[str, Dict[int, int]]
"""
{gene_type: {n_to_del: freq}}
"""

lh = get_logger(__name__)


class GenerationFailure(RuntimeError):
    ...


class Cdr3InsertionTable:
    _table: Dict[str, Dict[int, List[List[int]]]]
    """
    {chain: {num_to_insert: [[freq]]}}
    """
    _num_to_insert: Dict[str, Dict[int, int]]
    """
    {chain: {num_to_insert: freq}}
    """

    def __init__(self, cdr3_insertion_table_path: str) -> None:
        with get_reader(cdr3_insertion_table_path) as reader:
            self._table = json.load(reader)
        self._table["A"] = {
            int(k): v
            for k, v in self._table["A"].items()
        }
        self._table["B"] = {
            int(k): v
            for k, v in self._table["B"].items()
        }
        self._num_to_insert = {
            "A": {
                i: sum(itertools.chain(*self._table["A"][i]))
                for i in self._table["A"].keys()
            },
            "B": {
                i: sum(itertools.chain(*self._table["B"][i]))
                for i in self._table["B"].keys()
            }
        }

    def generate_cdr3(self, chain: str) -> TCRTranslationTableType:
        rets = []
        cdr3_length = random.choices(
            list(self._num_to_insert[chain].keys()),
            list(self._num_to_insert[chain].values())
        )[0]
        cdr3_pos_freq: List[List[int]] = self._table[chain][cdr3_length]
        for i in range(cdr3_length):
            aa = random.choices(
                self._table["AANames"],
                cdr3_pos_freq[i]
            )[0]
            aa_transl_table = TRANSL_TABLES[1]["AA"]
            rets.append((
                TRANSL_TABLES_NT[
                    random.choice([index for index in range(len(aa_transl_table)) if aa_transl_table[index] == aa])
                ],
                aa,
                aa,
                "="
            ))
        return rets


class TCell:
    _traj_name: str
    _trbj_name: str
    _trav_name: str
    _trbv_name: str
    _traj_tt: TCRTranslationTableType
    _trav_tt: TCRTranslationTableType
    _trbj_tt: TCRTranslationTableType
    _trbv_tt: TCRTranslationTableType
    _tra_cdr3_tt: TCRTranslationTableType
    _trb_cdr3_tt: TCRTranslationTableType
    _cell_uuid: str

    def __init__(
            self,
            cell_uuid: str,
            traj_name: str,
            trbj_name: str,
            trav_name: str,
            trbv_name: str,
            traj_tt: TCRTranslationTableType,
            trav_tt: TCRTranslationTableType,
            trbj_tt: TCRTranslationTableType,
            trbv_tt: TCRTranslationTableType,
            tra_cdr3_tt: TCRTranslationTableType,
            trb_cdr3_tt: TCRTranslationTableType
    ):
        self._cell_uuid = cell_uuid
        self._traj_name = traj_name
        self._trbj_name = trbj_name
        self._trav_name = trav_name
        self._trbv_name = trbv_name
        self._traj_tt = traj_tt
        self._trav_tt = trav_tt
        self._trbj_tt = trbj_tt
        self._trbv_tt = trbv_tt
        self._tra_cdr3_tt = tra_cdr3_tt
        self._trb_cdr3_tt = trb_cdr3_tt

    @classmethod
    def from_gene_names(
            cls,
            tcr_genelist: Dict[str, List[str]],
            cdr3_deletion_table: Cdr3DeletionTableType,
            cdr3_insertion_table: Cdr3InsertionTable,
            tcr_cache: Dict[str, TCRTranslationTableType]
    ):
        def choose_name(tcr_type: str) -> Tuple[str, TCRTranslationTableType]:
            while True:
                tcr_name = random.choice(tcr_genelist[tcr_type])
                if tcr_name in tcr_cache:
                    return tcr_name, tcr_cache[tcr_name]

        def clip_aa(tr_cdr3_tt, trv_tt, trj_tt):
            trv_tt_real_aa = "".join(list(zip(*trv_tt))[2])
            trj_tt_real_aa = "".join(list(zip(*trj_tt))[2])
            c_idx = trv_tt_real_aa[::-1].find("C")
            if c_idx == -1 or c_idx > 0.5 * len(trv_tt):
                c_idx = int(0.5 * len(trv_tt))
            f_idx = trj_tt_real_aa.find("F")
            if f_idx == -1 or f_idx > 0.5 * len(trj_tt):
                f_idx = int(0.5 * len(trj_tt))
            tr_cdr3_tt = trv_tt[- c_idx - 1: -1] + tr_cdr3_tt + trj_tt[0: f_idx + 1]
            trv_tt = trv_tt[0: - c_idx - 1]
            trj_tt = trj_tt[f_idx:]
            return tr_cdr3_tt, trv_tt, trj_tt

        (trbv_name, trbv_tt), (trbj_name, trbj_tt) = choose_name("trbv_names"), choose_name("trbj_names")
        (traj_name, traj_tt), (trav_name, trav_tt) = choose_name("traj_names"), choose_name("trav_names")

        chosen_deletion: Dict[str, int] = {
            k: random.choices(
                population=list(cdr3_deletion_table[k].keys()),
                weights=list(cdr3_deletion_table[k].values())
            )[0]
            for k in cdr3_deletion_table.keys()
        }

        def pop(tr_tt, gene_name, pos):
            while tr_tt[pos][-1] == "X":
                _ = tr_tt.pop(pos)
            for _ in range(chosen_deletion[gene_name]):
                _ = tr_tt.pop(pos)
            return tr_tt

        try:
            trav_tt = pop(trav_tt, "trav", -1)
            trbv_tt = pop(trbv_tt, "trbv", -1)
            traj_tt = pop(traj_tt, "traj", 0)
            trbj_tt = pop(trbj_tt, "trbj", 0)
        except (IndexError, RuntimeError) as e:
            raise GenerationFailure from e

        tra_cdr3_tt = cdr3_insertion_table.generate_cdr3("A")
        trb_cdr3_tt = cdr3_insertion_table.generate_cdr3("B")

        try:
            tra_cdr3_tt, trav_tt, traj_tt = clip_aa(tra_cdr3_tt, trav_tt, traj_tt)
            trb_cdr3_tt, trbv_tt, trbj_tt = clip_aa(trb_cdr3_tt, trbv_tt, trbj_tt)
        except (IndexError, RuntimeError) as e:
            raise GenerationFailure from e

        if len(trbv_tt) * len(trav_tt) * len(trbj_tt) * len(traj_tt) == 0:
            raise GenerationFailure

        return cls(
            cell_uuid=str(uuid.uuid4()),
            trav_tt=trav_tt,
            traj_tt=traj_tt,
            tra_cdr3_tt=tra_cdr3_tt,
            trbv_tt=trbv_tt,
            trbj_tt=trbj_tt,
            trb_cdr3_tt=trb_cdr3_tt,
            traj_name=traj_name,
            trav_name=trav_name,
            trbj_name=trbj_name,
            trbv_name=trbv_name
        )

    def to_fasta_record(self) -> str:
        return "\n".join((
            f">{self._cell_uuid}:A",
            self.alpha_nt,
            f">{self._cell_uuid}:B",
            self.beta_nt
        ))

    @property
    def cell_uuid(self):
        return self._cell_uuid

    @property
    def alpha_names(self) -> Tuple[str, str]:
        return self._trav_name, self._traj_name

    @property
    def beta_names(self) -> Tuple[str, str]:
        return self._trbv_name, self._trbj_name

    @property
    def alpha_nt(self) -> str:
        return "".join(
            itertools.chain(
                list(zip(*self._traj_tt))[0],
                list(zip(*self._tra_cdr3_tt))[0],
                list(zip(*self._trav_tt))[0]
            )
        )

    @property
    def beta_nt(self) -> str:
        return "".join(
            itertools.chain(
                list(zip(*self._trbj_tt))[0],
                list(zip(*self._trb_cdr3_tt))[0],
                list(zip(*self._trbv_tt))[0]
            )
        )

    @property
    def cdr3_aa(self) -> Tuple[str, str]:
        return (
            "".join(list(zip(*self._tra_cdr3_tt))[1]),
            "".join(list(zip(*self._trb_cdr3_tt))[1])
        )


def rearrange_tcr(
        n_cells: int,
        tcr_cache_path: str,
        tcr_genelist_path: str,
        output_base_path: str,
        cdr3_deletion_table_path: str,
        cdr3_insertion_table_path: str
):
    n_failure = 0
    with get_reader(tcr_genelist_path) as reader:
        tcr_genelist = json.load(reader)
    with get_reader(cdr3_deletion_table_path) as reader:
        cdr3_deletion_table = json.load(reader)
    with get_reader(tcr_cache_path) as reader:
        tcr_cache: Dict[str, TCRTranslationTableType] = json.load(reader)
    cdr3_deletion_table = {
        k: {
            int(_k): _v
            for _k, _v in v.items()
        }
        for k, v in cdr3_deletion_table.items()
    }
    with get_writer(output_base_path + ".fa") as writer:
        with load_table_appender_class("TSVTableAppender")(
                filename=output_base_path + ".stats",
                header=[
                    "UUID", "TRAV", "TRAJ", "TRBV", "TRBJ", "ACDR3_AA", "BCDR3_AA"
                ],
                tac=TableAppenderConfig(buffer_size=1024)
        ) as appender:
            for _ in tqdm(range(n_cells), desc="Simulating..."):
                while True:
                    try:
                        cell = TCell.from_gene_names(
                            tcr_genelist=tcr_genelist,
                            cdr3_deletion_table=cdr3_deletion_table,
                            cdr3_insertion_table=Cdr3InsertionTable(cdr3_insertion_table_path),
                            tcr_cache=tcr_cache
                        )
                    except GenerationFailure:
                        n_failure += 1
                        continue
                    else:
                        break
                writer.write(
                    cell.to_fasta_record() + "\n"
                )
                appender.append([
                    cell.cell_uuid,
                    *cell.alpha_names,
                    *cell.beta_names,
                    *cell.cdr3_aa
                ])
    lh.info("Finished with %d failures", n_failure)


if __name__ == "__main__":
    basename = "."
    rm_rf(os.path.join(basename, "sim") + ".fa.d")
    rearrange_tcr(
        n_cells=100,
        output_base_path=os.path.join(basename, "sim"),
        tcr_genelist_path=os.path.join(basename, "tcr_genelist.json"),
        tcr_cache_path=os.path.join(basename, "tcr_cache.json"),
        cdr3_deletion_table_path=os.path.join(basename, "cdr3_deletion_table.json"),
        cdr3_insertion_table_path=os.path.join(basename, "cdr3_insertion_table.json")
    )
