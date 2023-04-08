#!/usr/bin/env bash
# axel https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 
# axel https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
# gunzip ./*.gz
# samtools faidx hg38.fa

# python -m yasim transcribe -f hg38.fa -g hg38.ncbiRefSeq.gtf -o hg38_trans.fa

function perform_housekeeping() {
    rm -rf "${1}".d && \
    touch "${1}".finished || return 1
}
function perform_pbsim2_simulation() {
    OUTPUT_BASENAME="${1}"_pbsim2_"${2}"
    LOG_FILE_NAME="yasim_${OUTPUT_BASENAME}.log" \
        python -m yasim pbsim2 \
        -e ../bin/pbsim2 \
        -m "${2}" \
        -F hg38_trans.fa.d \
        -d "${1}".tsv \
        -o "${OUTPUT_BASENAME}" \
        -j 40 || return
    perform_housekeeping "${OUTPUT_BASENAME}"
}
perform_pbsim2_simulation FINE2a_chr1.depth R103 &
perform_pbsim2_simulation FINE2b_chr1.depth R103 &
perform_pbsim2_simulation TesR7A_chr1.depth R103 &
perform_pbsim2_simulation TesR7B_chr1.depth R103 &
wait
