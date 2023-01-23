#!/usr/bin/env bash
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/genes/ce11.ncbiRefSeq.gtf.gz
axel https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz

if ! mamba env list | grep ^yasim-pbsim1 &>> /dev/null; then
    mamba create -y -n yasim-pbsim1 -c bioconda pbsim=1.0.3
fi
if ! mamba env list | grep ^yasim-pbsim2 &>> /dev/null; then
    mamba create -y -n yasim-pbsim2 -c bioconda pbsim2=2.0.1
fi
if ! mamba env list | grep ^yasim-badread &>> /dev/null; then
    mamba create -y -n yasim-badread -c bioconda badread=0.2.0 python-edlib
fi

PYTHONPATH=../../src python -m yasim generate_as_events -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq_as.gtf
PYTHONPATH=../../src python -m yasim generate_depth -g ce11.ncbiRefSeq_as.gtf -o ce11_depth.tsv -d 20
PYTHONPATH=../../src python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq_as.gtf -o ce11_trans_as.fa
python -m labw_utils.bioutils get_gtf_statistics -g ce11.ncbiRefSeq_as.gtf -o ce11.ncbiRefSeq_as.gtf
PYTHONPATH=../../src python -m yasim pbsim -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/pbsim.sh" -j 40 -o ce11_pbsim_clr
PYTHONPATH=../../src python -m yasim pbsim2 -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/pbsim2.sh" -j 40 -o ce11_pbsim2_r94 -m R94
PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_nanopore2018 -m nanopore2018
PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_pacbio2016 -m pacbio2016
PYTHONPATH=../../src python -m yasim badread -d ce11_depth.tsv -F ce11_trans_as.fa.d -e "wrappers/badread.sh" -j 40 -o ce11_badread_nanopore2020 -m nanopore2020
python -m labw_utils.bioutils describe_fastq ce11_badread_nanopore2018.fq ce11_badread_nanopore2020.fq ce11_badread_pacbio2016.fq ce11_pbsim2_r94.fq ce11_pbsim_clr.fq

for fn in *.fq.stats.d; do
    Rscript R/plot_describe_fastq.R "${fn}"
done

PYTHONPATH=../../src python -m yasim transcribe -f ce11.fa -g ce11.ncbiRefSeq.gtf -o ce11_trans.fa
python -m labw_utils.bioutils get_gtf_statistics -g ce11.ncbiRefSeq.gtf -o ce11.ncbiRefSeq.gtf

STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ce11_STAR_index --genomeFastaFiles ce11.fa --sjdbGTFfile ce11.ncbiRefSeq.gtf
# !!!!! WARNING: --genomeSAindexNbases 14 is too large for the genome size=100286401, which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases 12
mkdir -p ce11_trans_bwa_index/
bwa index ce11_trans.fa -p ce11_trans_bwa_index/ce11

for fn in *.fq; do
    minimap2 -x splice -a -t 50 ce11.fa "${fn}" > "${fn}".sam
    samtools sort "${fn}".sam -@ 50 -o "${fn}".bam
    samtools index "${fn}".bam
    rm "${fn}".sam
done

python -m labw_utils.bioutils describe_sam ce11_badread_nanopore2018.fq.bam ce11_badread_nanopore2020.fq.bam ce11_badread_pacbio2016.fq.bam ce11_pbsim2_r94.fq.bam ce11_pbsim_clr.fq.bam

for fn in *.fq; do
    minimap2 -x splice -a -t 50 ce11_trans.fa "${fn}" > "${fn/.fq/_trans.fq}".sam
    samtools sort "${fn/.fq/_trans.fq}".sam -@ 50 -o "${fn/.fq/_trans.fq}".bam
    samtools index "${fn/.fq/_trans.fq}".bam
    rm "${fn/.fq/_trans.fq}".sam
done

python -m labw_utils.bioutils describe_sam ce11_badread_nanopore2018_trans.fq.bam ce11_badread_nanopore2020_trans.fq.bam ce11_badread_pacbio2016_trans.fq.bam ce11_pbsim2_r94_trans.fq.bam ce11_pbsim_clr_trans.fq.bam

for fn in *.fq.bam.stats.d; do
    Rscript R/plot_describe_sam.R "${fn}"
done

for fn in *trans.fq.bam.stats.d; do
    Rscript R/transform_depth_results.R "${fn}"/pileup_stat.tsv.gz "${fn}"/pileup_stat.merged.tsv
done

for fn in ce11_badread_nanopore2018.fq.bam ce11_badread_nanopore2020.fq.bam ce11_badread_pacbio2016.fq.bam ce11_pbsim2_r94.fq.bam ce11_pbsim_clr.fq.bam; do
     featureCounts -O -M --primary --ignoreDup -a ce11.ncbiRefSeq.gtf -g transcript_id -L -o "${fn}".fc.tsv "${fn}"
     stringtie -L -G ce11.ncbiRefSeq.gtf -o "${fn}".stringtie.gtf -p 50 "${fn}"
    python -m labw_utils.bioutils get_gtf_statistics -g "${fn}".stringtie.gtf -o "${fn}".stringtie.gtf
done
