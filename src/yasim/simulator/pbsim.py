import os
from typing import List

from commonutils import ioctl
from yasim.simulator import Simulator

FILE_DIR = os.path.dirname(__file__)

#     pbsim --prefix "${SIM_DIR}"/pbsim/"${2}"/sd \
#     --data-type CLR  "${1}" \
#     --model_qc /usr/share/pbsim/models/model_qc_clr \
#     --length-mean 3080 \
#     --length-sd 2211 \
#     --length-min 50 \
#     --length-max 50000 \
#     --accuracy-mean 1.0 \
#     --accuracy-sd 0.0 \
#     --accuracy-min 1.0 \
#     --difference-ratio 47:38:15 \
#     "${SIM_DIR}"/hg38_shuffled_filtered_"${2}".fa

class _SimulatorPbsimBase(Simulator):
    tmp_prefix: str

    def assemble_cmd(self) -> List[str]:
        cmd = [
            "pbsim",
            "--prefix",
            self.tmp_prefix,
            "--depth",
            str(self.depth),
            self.input_fasta
        ]
        return cmd

    def move_file_after_finish(self):
        pass # TODO

    def run(self) -> None:
        self.run_simulator_as_process("pbsim")

class SimulatePbsimCLR(_SimulatorPbsimBase):
    def assemble_cmd(self) -> List[str]:
        cmd = [
            "pbsim",
            "--prefix",
            self.tmp_prefix,
            "--depth",
            str(self.depth),
            "--data-type",
            "CLR",
            "--model_qc",
            os.path.join(FILE_DIR, "pbsim_dist", "model_qc_clr"),
            self.input_fasta
        ]
        return cmd


class SimulatePbsimCCS(_SimulatorPbsimBase):
    def assemble_cmd(self) -> List[str]:
        cmd = [
            "pbsim",
            "--prefix",
            self.tmp_prefix,
            "--depth",
            str(self.depth),
            "--data-type",
            "CCS",
            "--model_qc",
            os.path.join(FILE_DIR, "pbsim_dist", "model_qc_ccs"),
            self.input_fasta
        ]
        return cmd

