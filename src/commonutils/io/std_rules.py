import bz2
import gzip
import lzma

from commonutils.io.rules import BaseFileRule, FileRuleRing


def use_python_std_rules():
    FileRuleRing.register(
        BaseFileRule.from_extension(".gz", ".GZ", opener=gzip.open, rule_name="gzip_suffix_rule")
    )
    FileRuleRing.register(
        BaseFileRule.from_extension(".xz", ".lzma", opener=lzma.open, rule_name="lzma_suffix_rule")
    )
    FileRuleRing.register(
        BaseFileRule.from_extension(".bz2", opener=bz2.open, rule_name="bz2_suffix_rule")
    )
