# Batched C. Elegans Data Samples

## Introduction

Here lists datasets used for analysis of GEP, etc. The datasets are categorized by its sequencers (aka., instruments) and organized into following form:

- [ENA ACCESSION](about:blank)
  - Description about this project.
  - Date last updated of metadata by YASIM.
  - SINGLE or PAIRED or others.
  - Dataset being chosen (sample name used):
    - ENA ACCESSION (sample name), SIZE
  - Reference.

Some sample may have `**MAY NOT BE USABLE -- REQUIRE REVIEW**`. This means whether to use them is not sure.

## Illumina

### HiSeq 1500

- [PRJNA726525](https://www.ebi.ac.uk/ena/browser/view/PRJNA726525)
  - We report that CGP37157 improved the evolution with age of the sarcomeric regular structure, delaying development of sarcopenia in C. elegans body wall muscle. Similarly, CGP37157 favoured the maintenance of a regular mitochondrial structure during aging. CGP37157 induced a four-fold increase in the expression of ncx-6, one of the C. elegans mitochondrial Na+/Ca2+ exchangers. Overall design: Examination of 3 replicates of CGP37157 treated C.elegans and 3 replicates of control C.elegans.
  - Updated 2022-04-14.
  - paired
  - data (WT):
    - SRR14369911 (WT C. elegans control rep.1 \[Control_rep1_02_12\]), 1.1 GB 1.2 GB
    - SRR14369912 (WT C. elegans control rep.2 \[Control_rep2_04_12\]), 1.1 GB 1.2 GB
    - SRR14369913 (WT C. elegans control rep.3 \[Control_rep3_17_12\]), 1.1 GB 1.2 GB
  - The Mitochondrial Na<sup>+</sup>/Ca<sup>2+</sup> Exchanger Inhibitor CGP37157 Preserves Muscle Structure and Function to Increase Lifespan and Healthspan in Caenorhabditis elegans. García-Casas P, Alvarez-Illera P, Gómez-Orte E, Cabello J, Fonteriz RI, Montero M, Alvarez J. Departamento de Bioquímica y Biología Molecular y Fisiología, Facultad de Medicina, Unidad de Excelencia Instituto de Biología y Genética Molecular (IBGM), Universidad de Valladolid and CSIC, Valladolid, Spain. Front Pharmacol 12695687 (2021 ) <https://doi.org/10.3389/fphar.2021.695687>

### HiSeq 2000

- [PRJNA321853](https://www.ebi.ac.uk/ena/browser/view/PRJNA321853)
  - To understand the function and regulation of the C. elegans heat shock factor (HSF-1) in larval development, we have used ChIP-seq to analyze the occupancy of HSF1 and RNA Pol II in L2 larvae and young adult (YA) animals grown at 20°C or upon heat shock at 34°C for 30 min. In addition, we have used RNA-seq to analyze the transcriptomes of wild type (N2), hsf-1(ok600) mutants and hsf-1(ok600); rmSi1\[hsf-1::gfp\] L2 larvae grown at 20°C and characterized the gene expression change by heat shock in wild type (N2), hsf-1(sy441) and hsf-1(sy441);rmSi1\[hsf-1::gfp\] animals at L2 stage. Overall design: Experiment type: RNA-seq. Biological Source: strain: N2, hsf-1(sy441), hsf-1(sy441);rmSi1\[hsf-1::gfp\]; developmental dtage: L2 Larva. Experimental Factors: temperature: 20 degree celsius, 34 degree celsius.
  - Updated 2022-04-14.
  - single
  - data (WT):
    - SRR3535753 (RNA-seq of L2 larvae at 20°C, replicate 1), 906.4 MB
    - SRR3535754 (RNA-seq of L2 larvae at 20°C, replicate 2), 1.1 GB
    - SRR3535755 (RNA-seq of L2 larvae at 20°C, replicate 3), 824.3 MB
  - E2F coregulates an essential HSF developmental program that is distinct from the heat-shock response. Li J, Chauve L, Phelps G, Brielmann RM, Morimoto RI. Department of Molecular Biosciences, Rice Institute for Biomedical Research, Northwestern University, Evanston, Illinois 60208, USA. Genes Dev 30(18): 2062-2075 (2016 Sep) <https://doi.org/10.1101/gad.283317.116>

### HiSeq 2500

- [PRJNA412861](https://www.ebi.ac.uk/ena/browser/view/PRJNA412861)
  - RNA-seq libraries were generated for L4 wildtype worms to analyze RNA abundance at this point in development. Overall design: Three independent replicates of wildtype C. elegans were collected at the L4 stage for RNA.
  - Updated 2022-04-14.
  - single
  - data:
    - SRR6125151 (WT rep1 (RNA-Seq)), 1.9 GB
    - SRR6125152 (WT rep2 (RNA-Seq)), 1.7 GB
    - SRR6125153 (WT rep3 (RNA-Seq)), 1.7 GB
  - Short poly(A) tails are a conserved feature of highly expressed genes. Lima SA, Chipman LB, Nicholson AL, Chen YH, Yee BA, Yeo GW, Coller J, Pasquinelli AE. Division of Biology, University of California, San Diego, La Jolla, California, USA. Nat Struct Mol Biol 24(12): 1057-1063 (2017 Dec) <https://doi.org/10.1038/nsmb.3499>
- [PRJNA715371](https://www.ebi.ac.uk/ena/browser/view/PRJNA715371)
  - Bulk RNA seq of FACS isolated C. elegans neurons, with pan-neuronal reference, and sorted viable cell reference samples. Collected for comparison to single cell sequencing data Overall design: Cell type specific C. elegans neuronal samples were obtained via FACS for comparison to each other, and to pan-neuronal and all cell references. Sample names describe the cell type targeted (ASG, AVE, AWA, etc…), Rab3 indicates samples targeted with a Pan-Neuronal Marker for direct comparison to each individual cell type. Ref samples were obtained by dissociating and FACS isolation using only stains to sepatate living and dead cells and thus represents an "all cell" reference.
  - ** TOO SPECIFIC AND TOO SMALL **
  - Updated 2022-04-14.
  - paired
  - data (ASGr):
    - SRR13995310 (ASGr13), 486.6 MB 530.7 MB
    - SRR13995311 (ASGr14), 392.8 MB 434.2 MB
    - SRR13995312 (ASGr15), 414.0 MB 455.7 MB
    - SRR13995313 (ASGr16), 401.0 MB 453.0 MB
  - A head-to-head comparison of ribodepletion and polyA selection approaches for C. elegans low input RNA-sequencing libraries. Barrett A, McWhirter R, Taylor SR, Weinreb A, Miller DM, Hammarlund M. Department of Genetics, Yale University School of Medicine, New Haven, CT, 06510, USA. G3 (Bethesda) (2021 Apr) <https://doi.org/10.1093/g3journal/jkab121>
- [PRJNA549257](https://www.ebi.ac.uk/ena/browser/view/PRJNA549257)
  - Three independent replicates of wildtype C. elegans were collected at the L4 stage for RNA. Heat Shock was performed at 35C for 4hrs. Overall design: RNA-seq libraries were generated for L4 wildtype worms to analyze RNA expression in control vs HS worms.
  - Updated 2022-04-14.
  - paired
  - data (WT)
    - SRR9310678 (CTRL Rep 1 RNA-seq), 2.2 GB 2.3 GB
    - SRR9310679 (CTRL Rep 2 RNA-seq), 2.0 GB 2.1 GB
    - SRR9310680 (CTRL Rep 3 RNA-seq), 2.1 GB 2.1 GB
  - Remodeling of the Caenorhabditis elegans non-coding RNA transcriptome by heat shock. Schreiner WP, Pagliuso DC, Garrigues JM, Chen JS, Aalto AP, Pasquinelli AE. Division of Biology, University of California, San Diego, La Jolla, CA 92093-0349, USA. Nucleic Acids Res 47(18): 9829-9841 (2019 Oct) <https://doi.org/10.1093/nar/gkz693>

### HiSeq 3000

NO DATA!

### HiSeq 4000

- [PRJNA684142](https://www.ebi.ac.uk/ena/browser/view/PRJNA684142)
  - prp-40 regulates alternative splicing. In particular, microexons are highly dependent on prp-40 for inclusion.
  - **TOO BIG**
  - Updated 2022-04-14.
  - paired
  - data (WT):
    - SRR13238604 (wild type c), 10.9 GB 11.7 GB
    - SRR13238605 (wild type b), 5.3 GB 5.6 GB
    - SRR13238606 (wild type a), 3.2 GB 3.4 GB
  - NO REFERENCE
- [PRJNA803492](https://www.ebi.ac.uk/ena/browser/view/PRJNA803492)
  - We report the gene expression in L1 arrest (48 hours) of wild-type sleeping worms (N2) and aptf-1(gk794) mutant sleepless worms (HBR227) Overall design: There are 8 samples in total. 4 biological replicates of N2 and 4 biological replicates of HBR227
  - Updated 2022-04-14.
  - paired
  - data (WT):
    - SRR17886746 (N2-D), 2.3 GB 2.3 GB
    - SRR17886747 (N2-C), 2.0 GB 2.1 GB
    - SRR17886748 (N2-B), 2.4 GB 2.5 GB
    - SRR17886749 (N2-A), 2.0 GB 1.9 GB
  - Sleep neuron depolarization promotes protective gene expression changes and FOXO activation. Koutsoumparis A, Welp LM, Wulf A, Urlaub H, Meierhofer D, Börno S, Timmermann B, Busack I, Bringmann H. Chair of Cellular Circuits and Systems, Biotechnology Center (BIOTEC), Center for Molecular and Cellular Bioengineering (CMCB), Technical University Dresden, 01307 Dresden, Germany. Curr Biol 32(10): 2248-2262.e9 (2022 May) <https://doi.org/10.1016/j.cub.2022.04.012>



### HiSeq X

NO DATA!

### NextSeq 500

- [PRJNA448043](https://www.ebi.ac.uk/ena/browser/view/PRJNA448043)
  - The goal of RNA-seq is to identify the differential expressed genes in the wild-type worm and nmad-1 mutant worm at 20°C and 25°C respectively. 100 pairs of germlines were extracted and collected from worm tips. Two biological replicates were assigned for each group and the RNA profiles were generated by deep sequencing, using NEXTseq 500 Illumina. We mapped about 20 - 30 million reads per sample to C. elegans genome (build ce10) with bowtie2 workflow. Differential expression genes were identified through EBSeq in RSEM pipeline with 0.05 p value and 1.2 fold change. Our results showed elevated number of misregulated genes in the nmad-1 mutant groups at 25°C. Overall design: Germline mRNA profiles of wilde-type (WT) worm and nmad-1 mutant at 20°C and 25°C were generated by deep sequencing, in duplicate, using Nextseq500
  - Updated 2022-04-14.
  - paired
  - data (WT)
    - SRR6914611 (n2 rep1 20°C), 672.3 MB 710.6 MB
    - SRR6914612 (n2 rep2 20°C), 1.1 GB 1.1 GB
  - The demethylase NMAD-1 regulates DNA replication and repair in the Caenorhabditis elegans germline. Wang SY, Mao H, Shibuya H, Uzawa S, O'Brown ZK, Wesenberg S, Shin N, Saito TT, Gao J, Meyer BJ, Colaiácovo MP, Greer EL. Division of Newborn Medicine, Children's Hospital Boston, Boston, Massachusetts, United States of America. PLoS Genet 15(7): e1008252 (2019 Jul) <https://doi.org/10.1371/journal.pgen.1008252>
- [PRJNA488660](https://www.ebi.ac.uk/ena/browser/view/PRJNA488660)
  - RNA-seq of Wild Type (N2), pmk-1 or atf-7 mutant animals exposed to either non-pathogenic E. coli OP50 or pathogenic P. aeruginosa PA14 Overall design: mRNA profiles were generated using 3 replicates (>1,000 animals each) of each condition were prepared and sequenced, except for atf-7(qd22qd130) on PA14 which had only 2 replicates. Sequenced on Illumina NextSeq 500
  - Updated 2022-04-14.
  - paired
  - data (WT)
    - SRR7771392 (N2 OP50 1), 878.8 MB 914.8 MB
    - SRR7771394 (N2 OP50 2), 724.1 MB 788.0 MB
    - SRR7771396 (N2 OP50 3), 715.0 MB 750.0 MB
  - Global transcriptional regulation of innate immunity by ATF-7 in C. elegans. Fletcher M, Tillman EJ, Butty VL, Levine SS, Kim DH. Department of Biology, Massachusetts Institute of Technology, Cambridge, Massachusetts, United States of America. PLoS Genet 15(2): e1007830 (2019 Feb) <https://doi.org/10.1371/journal.pgen.1007830>

### NextSeq 550

<!-- - [PRJNA497368](https://www.ebi.ac.uk/ena/browser/view/PRJNA497368)
  - Small interfering RNAs (siRNAs) and their partner Argonaute proteins regulate the expression of target RNAs. When sperm and egg meet upon fertilization, a diverse set of proteins and RNA, including siRNA-Argonaute complexes, is passed on to the developing progeny. Thus, these two players are important to initiate specific gene expression programs in the next generation. The nematode Caenorhabditis elegans expresses several classes of siRNAs. 26G-RNAs are a particular class of siRNAs that are divided into two subpopulations: one expressed in the spermatogenic gonad and another expressed in oocytes and embryos. In this work, we describe the dynamics whereby oogenic 26G-RNAs setup gene silencing in the next generation. We also show several ways that spermatogenic 26G-RNAs and their partner Argonautes, ALG-3 and ALG-4, use to regulate their targets. Finally, we show that ALG-3 and ALG-4 are fine-tuning their own expression, a rare role of Argonaute proteins. Altogether, we provide new insights into how siRNAs and Argonautes are regulating gene expression.
  - 
  - Updated 2022-04-14.
  - single
  - data (WT)
    - SRR8081024 (Wild-type rep2), 488.9 MB
    - SRR8081025 (Wild-type rep1), 450.4 MB
    - SRR8081026 (Wild-type rep4), 486.2 MB
    - SRR8081027 (Wild-type rep3), 564.4 MB
  - Genome sequence of the nematode C. elegans: a platform for investigating biology. C. elegans Sequencing Consortium. Science 282(5396): 2012-2018 (1998 Dec) <https://doi.org/10.1126/science.282.5396.2012> -->

### NextSeq 1000

NO DATA!

### NextSeq 2000

NO DATA!

### NovaSeq 6000

- [PRJNA640412](https://www.ebi.ac.uk/ena/browser/view/PRJNA640412)
  - We analyzed RNA-seq data from C. elegans on six bacterial diets, 3 E. coli diets (OP50, HT115, HB101) found in the laboratory setting and 3 found in C. elegans natural environment (Methylobacterium, Xanthomonas, Sphingomonas). We compared gene expression between L4 worms on the OP50 diet to the other 5 bacterial diets. We found that diet influences many gene expression changes, leading to genes being both up regulated and down regulated in a diet-dependent manner. Overall design: mRNA profiles of L4 wildtype C. elegans on 6 different bacterial diets
  - Updated 2022-04-14.
  - paired
  - data (HB101, rep 1-3):
    - SRR12047107 (WT L4 worms fed HB101_rep1), 1.2 GB 1.3 GB
    - SRR12047108 (WT L4 worms fed HB101_rep2), 1.3 GB 1.4 GB
    - SRR12047108 (WT L4 worms fed HB101_rep3), 2.0 GB 2.0 GB
  - Bacterial diets differentially alter lifespan and healthspan trajectories in C. elegans. Stuhr NL, Curran SP. Leonard Davis School of Gerontology, University of Southern California, 3715 McClintock Ave, Los Angeles, CA, 90089, USA. Commun Biol 3(1): 653 (2020 Nov) <https://doi.org/10.1038/s42003-020-01379-1>

## AB SOLiD

NO DATA!

## Ion Torrent

NO DATA!

## Roche LS454

### 454 GS FLX Titanium

- [SRP007195](https://www.ebi.ac.uk/ena/browser/view/SRP007195)
  - A differential sequencing-based analysis of the C. elegans noncoding transcriptome
  - Updated 2022-04-14.
  - single
  - data
    - SRR278717 (C.elegans mixed mixed-stage worms), 33.7 MB
    - SRR278718 (C.elegans mixed mixed-stage worms), 35.6 MB
  - NO REFERENCE

## Helicos

### Helicos HeliScope

NO DATA!

## Oxford Nanopore

### MinION

- [PRJNA591184](https://www.ebi.ac.uk/ena/browser/view/PRJNA591184)
  - Despite highly conserved chromatin states and cis-regulatory elements, studies of metazoan genomes reveal that gene organization and the strategies to control mRNA expression can vary widely among animal species. C. elegans gene regulation is often assumed to be similar to that of other model organisms, yet evidence suggests the existence of distinct molecular mechanisms to pattern the developmental transcriptome, including extensive post-transcriptional RNA control pathways, widespread splice leader (SL) trans-splicing of pre-mRNAs, and the organization of genes into operons. Here, we performed ChIP-seq for histone modifications in highly synchronized embryos cohorts representing three major developmental stages, with the goal of better characterizing whether the dynamic changes in embryonic mRNA expression are accompanied by changes to the chromatin state. We were surprised to find that thousands of promoters are persistently marked by active histone modifications, despite a fundamental restructuring of the transcriptome. We employed global run-on sequencing using a long-read nanopore format to map nascent RNA transcription across embryogenesis, finding that the invariant open chromatin regions are persistently transcribed by Pol II at all stages of embryo development, even though the mature mRNA is not produced. By annotating our nascent RNA sequencing reads into directional transcription units, we find extensive evidence of polycistronic RNA transcription genome-wide, suggesting that nearby genes in C. elegans are linked by shared transcriptional regulatory mechanisms. We present data indicating that the sharing of cis-regulatory sequences has constrained C. elegans gene positioning and likely explains the remarkable retention of syntenic gene pairs over long evolutionary timescales. Overall design: Nascent Pol II transcription was analyzed by sequencing replicate nascent RNA libraries for gastrula stage and late stage C. elegans embryos; high-quality fastq reads from replicate Gro-seq libraries were combined into a single library representing gastrula or late embryo nascent transcriptomes
  - Updated 2022-04-14.
  - single
  - data
    - SRR10517637 (gastrula_GroSeq), 3.1 GB
    - SRR10517638 (late_emb_GroSeq), 5.0 GB
  - Genome sequence of the nematode C. elegans: a platform for investigating biology. C. elegans Sequencing Consortium. Science 282(5396): 2012-2018 (1998 Dec) <https://doi.org/10.1126/science.282.5396.2012>

- [PRJNA533634](https://www.ebi.ac.uk/ena/browser/view/PRJNA533634)
  - High throughput RNA sequencing (RNA-seq) using cDNA has played a key role in delineating transcriptome complexity, including alternative transcription initiation, splicing, polyadenylation and base modification. However, the reads derived from current RNA-seq technologies are usually short and deprived of information on modification during reverse transcription, compromising their potential in defining transcriptome complexity. Here we applied a direct RNA sequencing method with ultra-long reads from Oxford Nanopore Technologies (ONT) to study the transcriptome complexity in C. elegans. We sequenced native poly-A tailed mRNAs by generating approximately six million reads from embryos, L1 larvae and young adult animals, with average read lengths ranging from 900 to 1,100 bps across stages. Around half of the reads represent full-length transcripts, judged by the presence of a splicing-leader or their full coverage of an existing transcript. To take advantage of the full-length transcripts in defining transcriptome complexity, we devised a novel algorithm to predict novel isoforms or group them with exiting isoforms using their mapping tracks rather than the existing intron/exon structures, which allowed us to identify roughly 57,000 novel isoforms and recover at least 26,000 out of the 33,500 existing isoforms. Intriguingly, stage-specific expression at the level of gene and isoform demonstrates little correlation. Finally, we observed an elevated level of modification in all bases in the coding region relative to the UTR. Taken together, the ONT long reads are expected to deliver new insights into RNA processing and modification and their underlying biology. Overall design: Animals of different stages, i.e., embryo, L1 larva and young adult, were collected and total RNAs were extracted using TRIzol (Invitrogen) following the manufacturer’s instructions. Approximately 100 µg total RNAs were extracted for each sample. Around 900 ng of poly-A tailed mRNAs was purified using Dynabeads™ mRNA Purification Kit (Invitrogen) based on the user’s manual for each library preparation. Nanopore sequencing libraries were constructed using Direct RNA sequencing kit (cat# SQK-RNA001). The libraries were loaded onto Nanopore R9.4.1 flow cell (cat# FLO-MIN106) and sequenced on MinION acquired from Oxford Nanopore Technologies. The software used for sequencing was MINKNOW 2.1 with base-caller, Albacore (v2.0.1).
  - Updated 2022-04-14.
  - single
  - data:
    - SRR8929004 (Mix-stage embryo), 1.7 GB
    - SRR8929005 (L1 larva), 1.5 GB
    - SRR8929006 (Young adult), 2.0 GB
  - Genome sequence of the nematode C. elegans: a platform for investigating biology. C. elegans Sequencing Consortium. Science 282(5396): 2012-2018 (1998 Dec) <https://doi.org/10.1126/science.282.5396.2012>

### PromethION

NO DATA!

### GridION

NO DATA!

## PacBio

### PacBio Sequel

- [PRJNA764925](https://www.ebi.ac.uk/ena/browser/view/PRJNA764925)
  - To update the genome of a C. elegans Hawaiian strain, CB4856, and to identify presence-absence variants (PAVs) between two C. elegans strains, PD1074 and CB4856.
  - Updated 2022-04-14.
  - single
  - data:
    - SRR15993150 (CB4856_mixed_RNAseq_PacBio_Iso-Seq), 7.5 GB
    - SRR15993151 (PD1074_mixed_RNAseq_PacBio_Iso-Seq), 7.6 GB
  - Intraspecific <i>de novo</i> gene birth revealed by presence-absence variant genes in <i>Caenorhabditis elegans</i>. Lee BY, Kim J, Lee J. Research Institute of Basic Sciences, Seoul National University, Seoul 08826, Korea. NAR Genom Bioinform 4(2): lqac031 (2022 Jun) <https://doi.org/10.1093/nargab/lqac031>

### PacBio Sequel II

NO DATA!

## GenapSys

NO DATA!
