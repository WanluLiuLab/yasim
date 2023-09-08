import json
from labw_utils.commonutils.lwio.safe_io import get_reader
from labw_utils.typing_importer import List
from yasim.helper.translation_instruction import SimpleExon, SimpleTE, TranslationInstruction
from labw_utils.commonutils.appender import TableAppenderConfig, load_table_appender_class

import rmsk_parser

RMSK_OUT_GFF_PATH = "aln/ce11_denovo_test_rmsk_prep.rmsk.d/ce11_denovo_test_rmsk_prep.fa.out.gff"
RMSK_CONVERT_JSON_PATH = "sim/ce11_denovo_test_rmsk_prep.json"
GROUND_TRUTH_TI_PATH = "sim/ce11_denovo_test.json"
OUT_TSV_PATH = "aln/ce11_denovo_test_rmsk_prep.comp"


if __name__ == "__main__":
    with get_reader(RMSK_CONVERT_JSON_PATH, is_binary=False) as r:
        rmsk_conv: List[str] = json.load(r)
    with load_table_appender_class("TSVTableAppender")(
        filename=OUT_TSV_PATH,
        header=[
            "EVENT_ID",
            "ISOFORM_ID",
            "EVENT_TYPE",  # overlap_correct, unidentified, false_positive, overlap_incorrect
            "IDENTIFIED_TYPE",
            "CORRECT_TYPE",
            "OVERLAP_DICE",
        ],
        tac=TableAppenderConfig(),
    ) as appender:
        ti: TranslationInstruction = TranslationInstruction.from_json(GROUND_TRUTH_TI_PATH)
        all_te_positions = {}
        for ground_truth_transcript_name, ground_truth_transcript in ti.transcripts.items():
            te_positions = []
            curr_offset = 0
            for simple_segment in ground_truth_transcript.l:
                if isinstance(simple_segment, SimpleExon):
                    curr_offset += len(simple_segment.seq)
                elif isinstance(simple_segment, SimpleTE):
                    te_positions.append(
                        (curr_offset, curr_offset + len(simple_segment.seq) + 1, simple_segment.src_te_name)
                    )
                    curr_offset += len(simple_segment.seq)
            all_te_positions[ground_truth_transcript_name] = te_positions

        with rmsk_parser.RMSKGffIterator(RMSK_OUT_GFF_PATH, True) as rmski:
            last_feature_id = ""
            for feature in rmski:
                if (
                    feature.end0b - feature.start0b + 1 < 20
                    or "(" in feature.attribute_get("repeat_name")
                    or feature.attribute_get("repeat_name").endswith("-rich")
                ):
                    continue
                te_id_to_remove = -1
                current_feature_id = rmsk_conv[int(feature.seqname)]
                ground_truth_te_positions = all_te_positions[current_feature_id]
                for te_id, te in enumerate(ground_truth_te_positions):
                    if max(feature.start0b, te[0]) < min(feature.end0b, te[1]):
                        dice = (
                            (min(feature.end0b, te[1]) - max(feature.start0b, te[0]) + 1)
                            * 2
                            / (feature.end0b - feature.start0b + 1 + te[1] - te[0] + 1)
                        )
                        appender.append(
                            [
                                0,
                                feature.seqname,
                                "overlap_correct"
                                if feature.attribute_get("repeat_name") == te[2]
                                else "overlap_incorrect",
                                feature.attribute_get("repeat_name"),
                                te[2],
                                dice,
                            ]
                        )
                        if dice > 0.9:
                            te_id_to_remove = te_id
                        break
                else:
                    appender.append([0, feature.seqname, "false_positive", feature.attribute_get("repeat_name"), "", 0])
                if te_id_to_remove != -1:
                    ground_truth_te_positions.pop(te_id_to_remove)
                last_feature_id = current_feature_id
        for ground_truth_te_positions in all_te_positions.values():
            for te_id, te in enumerate(ground_truth_te_positions):
                appender.append([0, 0, "unidentified", "", te[2], 0])
