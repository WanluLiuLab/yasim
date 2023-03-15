#!/usr/bin/env bash
set -uex
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz
gunzip ./*.gz

grep -i '^chrI\s' < ce11.ncbiRefSeq.gtf > ce11.ncbiRefSeq.chr1.gtf
head ce11.fa -n "$(($(cat -n ce11.fa | grep '>' | head -n 2 | tail -n 1 | cut -f 1)-1))" >  ce11.chr1.fa

rm -f ce11.ncbiRefSeq.gtf ce11.fa

if ! mamba env list | grep ^yasim-pbsim1 &>> /dev/null; then
    mamba create -y -n yasim-pbsim1 -c bioconda pbsim=1.0.3
fi
if ! mamba env list | grep ^yasim-pbsim2 &>> /dev/null; then
    mamba create -y -n yasim-pbsim2 -c bioconda pbsim2=2.0.1
fi
if ! mamba env list | grep ^yasim-badread &>> /dev/null; then
    mamba create -y -n yasim-badread -c bioconda badread=0.2.0 python-edlib
fi

export PYTHONPATH=../../src

python -m yasim transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11_trans.chr1.fa
python -m labw_utils.bioutils describe_gtf ce11.ncbiRefSeq.chr1.gtf \
python -m yasim generate_as_events \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq.chr1.gtf \
    -o ce11.ncbiRefSeq_as.chr1.gtf \
    -c 2
python -m yasim transcribe \
    -f ce11.chr1.fa \
    -g ce11.ncbiRefSeq_as.chr1.gtf \
    -o ce11_trans_as.chr1.fa
python -m labw_utils.bioutils describe_gtf ce11.ncbiRefSeq_as.chr1.gtf

for dge in dge1 dge2; do
    python -m yasim generate_gene_depth \
        -g ce11.ncbiRefSeq_as.chr1.gtf \
        -o ce11_gene_depth_"${dge}".chr1.tsv \
        -d 20
    for diu in diu1 diu2; do
        python -m yasim generate_isoform_depth \
            -g ce11.ncbiRefSeq_as.chr1.gtf \
            -o ce11_depth_"${dge}"_"${diu}".chr1.tsv \
            -d ce11_gene_depth_"${dge}".chr1.tsv \
            --alpha 2
        python -m yasim generate_isoform_replicates \
            -d ce11_depth_"${dge}"_"${diu}".chr1.tsv \
            -n 2
    done
done

function perform_housekeeping() {
    gzip -9f "${1}".fq && \
    cat "${1}".d/*/*.maf | gzip -9f >"${1}".maf.gz && \
    rm -rf "${1}".d
}


function simulate(){
    python -m yasim pbsim \
        -d "${1}" \
        -F ce11_trans_as.chr1.fa.d \
        -e pbsim \
        -o ce11_pbsim_RS_clr_"${2}"
    perform_housekeeping ce11_pbsim_RS_clr_"${2}"

    python -m yasim pbsim \
        -d "${1}" \
        -F ce11_trans_as.chr1.fa.d \
        -e pbsim \
        --ccs \
        -o ce11_pbsim_RS_ccs_"${2}"
    perform_housekeeping ce11_pbsim_ccs_"${2}"

    for pbsim2_models in R95 P6C4; do
        python -m yasim pbsim2 \
        -d "${1}" \
        -F ce11_trans_as.chr1.fa.d \
        -e pbsim2 \
        -o ce11_pbsim2_"${pbsim2_models}"_"${2}" \
        -m "${pbsim2_models}"
        perform_housekeeping ce11_pbsim2_"${pbsim2_models}"_clr_"${2}"
    done

    python \
        -m yasim pbsim3 \
        -d "${1}" \
        -F ce11_trans_as.chr1.fa.d \
        -e pbsim3 \
        -o ce11_pbsim3_SEQUEL_clr_"${2}" \
        -m SEQUEL \
        -M errhmm
    perform_housekeeping ce11_pbsim3_SEQUEL_clr_"${2}"

    python \
        -m yasim pbsim3 \
        -d "${1}" \
        -F ce11_trans_as.chr1.fa.d \
        -e pbsim3 \
        -o ce11_pbsim3_SEQUEL_ccs_"${2}" \
        --ccs_pass 10 \
        -m SEQUEL \
        -M errhmm
    perform_housekeeping ce11_pbsim3_SEQUEL_ccs_"${2}"
}

# PBSIM 1, 2 and 3 is lightweighted so can be run parallely
for dge in dge1 dge2; do
    for diu in diu1 diu2; do
        for rep in 0 1; do
            simulate \
            ce11_depth_"${dge}"_"${diu}".chr1.tsv."${rep}" \
            "${dge}"_"${diu}"_"${rep}" &
        done
    done
done
wait

mkdir -p lastdb
lastdb -P40 lastdb/ce11_trans.chr1 ce11_trans.chr1.fa

for fn in *.fq.gz; do
    minimap2 -x splice -a -t 50 ce11.chr1.fa "${fn}" > "${fn}".sam
    samtools sort "${fn}".sam -@ 50 -o "${fn}".bam
    samtools index "${fn}".bam
    rm -f "${fn}".sam
    minimap2 -a -t 50 ce11_trans.chr1.fa "${fn}" > "${fn/.fq/_trans.fq}".sam
    samtools sort "${fn/.fq/_trans.fq}".sam -@ 50 -o "${fn/.fq/_trans.fq}".bam
    samtools index "${fn/.fq/_trans.fq}".bam
    rm -f "${fn/.fq/_trans.fq}".sam
    last-train -Qfastx -P40 lastdb/ce11_trans.chr1 "${fn}" > "${fn}"_trans.train
    lastal -P40 -Qfastx -m100 -j7 \
    -p "${fn}"_trans.train \
    lastdb/ce11_trans.chr1 "${fn}" |\
    last-map-probs /dev/stdin | \
    pigz -9 > "${fn}"_trans.maf
done

for fn in *.fq.gz; do
    python -m labw_utils.bioutils describe_fastq "${fn}" &
done
wait

for fn in *.fq.gz.bam; do
    python -m labw_utils.bioutils describe_sam "${fn}" &
done
wait

find . | grep .fq.gz.bam$ | grep -v trans | while read -r fn; do
    stringtie -L \
        -G ce11.ncbiRefSeq.chr1.gtf \
        -o "${fn}".stringtie.gtf \
        -p 40 \
        "${fn}" &
done
wait

for fn in *.stringtie.gtf; do
    python -m labw_utils.bioutils describe_gtf "${fn}" &
done
wait


find ./*.stringtie.gtf > stringtie-mergelist.txt
stringtie --merge -G ce11.ncbiRefSeq.chr1.gtf -p 40 -o stringtie_merged.gtf stringtie-mergelist.txt
python -m labw_utils.bioutils describe_gtf stringtie_merged.gtf
python -m labw_utils.bioutils transcribe -g stringtie_merged.gtf -f ce11.chr1.fa -o ce11_trans_stringtie.fa

find . | grep .fq.gz.bam$ | grep -v trans | while read -r fn; do
    featureCounts -L -O -M --primary --ignoreDup \
        -a ce11.ncbiRefSeq.chr1.gtf \
        -g transcript_id \
        -o "${fn}".fc.tsv \
        "${fn}" &
    featureCounts -L -O -M --primary --ignoreDup \
        -a ce11.ncbiRefSeq.chr1.gtf \
        -g gene_id \
        -o "${fn}".fc.gene.tsv \
        "${fn}" &
    featureCounts -L -O \
        -a ce11.ncbiRefSeq_as.chr1.gtf \
        -g transcript_id \
        -o "${fn}".fc_gt.tsv \
        "${fn}" &
    featureCounts -L -O \
        -a ce11.ncbiRefSeq_as.chr1.gtf \
        -g gene_id \
        -o "${fn}".fc_gt.gene.tsv \
        "${fn}" &
    featureCounts -L -O \
        -a stringtie_merged.gtf \
        -g transcript_id \
        -o "${fn}".fc_stringtie.tsv \
        "${fn}" &
    featureCounts -L -O \
        -a stringtie_merged.gtf \
        -g gene_id \
        -o "${fn}".fc_stringtie.gene.tsv \
        "${fn}" &
done
wait

find . | grep .fq.gz.bam$ | grep trans | while read -r fn; do
    {
        printf "REFERENCE_NAME\tPOS\tNUM_READS\n" > "${fn}".depth.tsv
        samtools depth -aa "${fn}" >> "${fn}".depth.tsv
        Rscript R/transform_depth_results.R "${fn}".depth.tsv "${fn}".depth.mean.tsv
    } &
done
wait

printf "FILENAME\tINSERTION\tDELETION\tMATCH\tSUBSTITUTION\n" > all_last_mapq.tsv
for fn in *.maf.gz; do
    python -m yasim_scripts extract_read_length_from_maf_yasim "${fn}" > "${fn}".rlen.tsv &
    python -m yasim_scripts extract_quality_from_maf "${fn}" >> all_last_mapq.tsv
done
wait

Rscript ./plot_fastq.R
Rscript ./plot_nipg.R
Rscript ./plot_sam.R
Rscript ./plot_gep.R
Rscript ./plot_maf_read_completeness.R
Rscript ./plot_maf_error_rate.R
