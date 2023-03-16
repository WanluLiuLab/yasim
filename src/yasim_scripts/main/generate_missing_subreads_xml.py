import glob
import os.path
import time
import uuid
from typing import List

import jinja2

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def main(args: List[str]):
    if len(args) != 1:
        _lh.error("Should have exact 1 arguments: The directory of CCS FASTQ")
        return
    env = jinja2.Environment(loader=jinja2.PackageLoader('yasim_scripts.main', 'templates'))
    template = env.get_template('pbsim_xml_template.xml')
    for fn in tqdm(glob.glob(os.path.join(args[0], "*.tmp.d", "tmp.subreads.bam"))):
        dataset_fn = fn+".xml"
        with get_writer(dataset_fn) as writer:
            timestamp = time.localtime()
            writer.write(template.render(
                timestamp_file=time.strftime("%y-%m-%dT%H:%M:%S", timestamp),
                timestamp_simple=time.strftime("%y%m%d_%H%m%S", timestamp),
                bam_filepath=fn,
                file_uuid=str(uuid.uuid4())
            ))

