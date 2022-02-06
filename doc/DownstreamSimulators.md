# Downstream Simulators

## Next Generation Sequencing


### 454sim
- C++
- Roche-454
- <https://sourceforge.net/projects/bioinfo-454sim/>
- Lysholm, F., Andersson, B., & Persson, B. (2011). An efficient simulator of 454 data using configurable statistical models. BMC research notes, 4, 449. <https://doi.org/10.1186/1756-0500-4-449>
- Scopus Cited 24

### PIRS
- C++, Perl
- Illumina
- <ftp://ftp.genomics.org.cn/pub/pIRS/>
- Xuesong Hu, Jianying Yuan, Yujian Shi, Jianliang Lu, Binghang Liu, Zhenyu Li, Yanxiang Chen, Desheng Mu, Hao Zhang, Nan Li, Zhen Yue, Fan Bai, Heng Li, Wei Fan, pIRS: Profile-based Illumina pair-end reads simulator, Bioinformatics, Volume 28, Issue 11, 1 June 2012, Pages 1533–1535, <https://doi.org/10.1093/bioinformatics/bts187>
- Scopus Cited 115

### flowsim
- Haskell, GPL
- Roche-454
- <https://hackage.haskell.org/package/flowsim>
- Balzer, S., Malde, K., Lanzén, A., Sharma, A., & Jonassen, I. (2010). Characteristics of 454 pyrosequencing data--enabling realistic simulation with flowsim. Bioinformatics (Oxford, England), 26(18), i420–i425. <https://doi.org/10.1093/bioinformatics/btq365>
- Scopus Cited 103

### ART
- C++
- Roche-454, Illumina, ABI-SOLiD
- <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>
- Huang, W., Li, L., Myers, J. R., & Marth, G. T. (2012). ART: a next-generation sequencing read simulator. Bioinformatics (Oxford, England), 28(4), 593–594. <https://doi.org/10.1093/bioinformatics/btr708>
- Scopus Cited 679

### Grinder
- Perl
- \*
- <https://sourceforge.net/projects/biogrinder/>
- Angly, F. E., Willner, D., Rohwer, F., Hugenholtz, P., & Tyson, G. W. (2012). Grinder: a versatile amplicon and shotgun sequence simulator. Nucleic acids research, 40(12), e94. <https://doi.org/10.1093/nar/gks251>
- Scopus Cited 129

### GemSIM
- Python
- Roche-454, Illumina
- <https://sourceforge.net/projects/gemsim>
- McElroy, K. E., Luciani, F., & Thomas, T. (2012). GemSIM: general, error-model based simulator of next-generation sequencing data. BMC genomics, 13, 74. <https://doi.org/10.1186/1471-2164-13-74>
- Scopus Cited 115

### XS
- C, GPL V3
- Ion Torrent, Roche-454, Illumina, ABI-SOLiD
- <https://bioinformatics.ua.pt/software/xs/>
- Pratas, D., Pinho, A. J., & Rodrigues, J. M. (2014). XS: a FASTQ read simulator. BMC research notes, 7, 40. <https://doi.org/10.1186/1756-0500-7-40>
- Scopus Cited 14

### curesim
- Java
- Ion Torrent
- <http://www.pegase-biosciences.com/curesim-a-customized-read-simulator/>
- Caboche, S., Audebert, C., Lemoine, Y., & Hot, D. (2014). Comparison of mapping algorithms used in high-throughput sequencing: application to Ion Torrent data. BMC genomics, 15, 264. <https://doi.org/10.1186/1471-2164-15-264>
- Scopus Cited 57

### FASTQSim
- Python and Java
- \*
- <https://sourceforge.net/p/fastqsim>
- Shcherbina A. (2014). FASTQSim: platform-independent data characterization and in silico read generation for NGS datasets. BMC research notes, 7, 533. <https://doi.org/10.1186/1756-0500-7-533>
- Scopus Cited 22

### dwgsim
- C
- Ion Torrent, Illumina and ABI-SOLiD
- <https://github.com/nh13/DWGSIM>

### wgsim
- C
- ?
- <https://github.com/lh3/wgsim>

### SimSeq
- C
- Illumina
- <https://github.com/jstjohn/SimSeq>

## Third Generation Sequencing

### pbsim

#### metadata
- Metadata:
  - C
  - PacBio CLR, PacBio CCS
  - <https://salsa.debian.org/med-team/pbsim>
  - Yukiteru Ono, Kiyoshi Asai, Michiaki Hamada, PBSIM: PacBio reads simulator—toward accurate genome assembly, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 119–121, <https://doi.org/10.1093/bioinformatics/bts649>
  - Scopus Cited 149
- Highlights:
  - Can perform sampling-based simulation (in which both length and quality scores are sampled from a real read set) and model-based simulation.
- Article:
  - PacBio have CLR (long reads with high error rates) and CCS (contrary) reads.
  - Length: log-normal distributions
  - Accuracy: normal distribution in CLR, exponential in CCS
  - Quality scores: no position-specific error profile in CLR and CCS reads was found
    - Simulated using a frequency table computed from _E. Coli_.
  - Simulation
    - Substitution and insertion simulated from quality scores. Bias found but not clear, so not introduced.
    - Deletion is uniform.
    - Coverage is uniform to GC content.

### SiLiCo

- Metadata:
  - Python
  - PacBio, Oxford Nanopore
  - <https://github.com/ethanagb/SiLiCO>
  - Baker, E. A. G., Goodwin, S., McCombie, W. R., & Mendivil Ramos, O. (2016). SiLiCO: A Simulator of Long Read Sequencing in PacBio and Oxford Nanopore [Preprint]. Genomics. <https://doi.org/10.1101/076901>
  - Scopus Cited ?
- Highlights:
  - pass
- Article:
  - Method:
    - A model-fitting approach was implemented on exemplary data sets from each platform.
    - Read Length:
      - PacBio results follow an approximately log-normal distribution with a mean read length of 3kb
      - Several candidate distributions (Weibull, gamma (best fits), log-normal) were fitted to Oxford Nanopore data sets.
    - Can perform Monte Carlo simulation.


### LongISLND

- Metadata:
  - Java, Python
  - Theoritically \*, tested on PacBio P5 and P6 chemistries
  - <https://github.com/bioinform/longislnd>
  - Lau, B., Mohiyuddin, M., Mu, J. C., Fang, L. T., Bani Asadi, N., Dallett, C., & Lam, H. Y. (2016). LongISLND: in silico sequencing of lengthy and noisy datatypes. Bioinformatics (Oxford, England), 32(24), 3829–3832. <https://doi.org/10.1093/bioinformatics/btw602>
  - Scopus Cited 11
- Highlights:
  - Can simulate FASTQ, H5 and BAM format
- Article:
  - Review:
    - pbsim
      - generates only the FASTQ data format
      - without the multi-pass mechanism 
      - without additional per-base probability and kinetic data required by downstream analysis tools
  - Method:
    - Independent of the underlying sequencing mechanism.
    - By a non-parametric learn-and-simulate approach.

### SimLoRD

- Metadata:
  - Python
  - PacBio CCS
  - <https://bitbucket.org/genomeinformatics/simlord/src/master/>
  - Stöcker, B. K., Köster, J., & Rahmann, S. (2016). SimLoRD: Simulation of Long Read Data. Bioinformatics (Oxford, England), 32(17), 2704–2706. <https://doi.org/10.1093/bioinformatics/btw286>
  - Scopus Cited ?
- Highlights:
  - pass 
- Article:
  - Review:
    - pbsim
      - Outdated chemistry and cannot be completely re-configured.
      - Conditional read quality distribution does not match well existing data.
      - does not provide SAM-formatted alignments between the reference and the simulated reads.
    - FASTQSim
      - Unable to provide mapping information or alignments of simulated reads
      - Simulates reads rather slowly
      - The simulated length/quality distributions do not agree well with data
      - It is difficult to change parameters directly. 
  - Error rate of CCSs decreases with the number of passes.
  - Read length: a log-normal distribution for genomic data and an empirical distribution corresponding to library size selection for RNA-seq data.
  - If the reference contains Ns in the relevant part, those Ns are replaced randomly in the read.
  - Error probabilities added by number of passes and read length from a scaled chi-squared distribution.

### NPBSS

- Metadata:
  - Matlab
  - PacBio CLR, PacBio CCS
  - <https://github.com/NWPU-903PR/NPBSS_Octave> \& <https://github.com/NWPU-903PR/NPBSS_MATLAB>
  - Wei, Z. G., & Zhang, S. W. (2018). NPBSS: a new PacBio sequencing simulator for generating the continuous long reads with an empirical model. BMC bioinformatics, 19(1), 177. <https://doi.org/10.1186/s12859-018-2208-0>
  - Scopus Cited 15
- Highlights:
  - pass
- Article:
  - Review
    - SimLoRD is more convenient than PBSIM and FASTQSim for parameters setting to simulate PacBio CCS reads
    - pbsim:
      - First, the quality value (QV, also called the Phred quality score) at each position for a simulated read is randomly chosen, but we found that the proportions of different QVs in real PacBio reads are different. 
      - Second, we also observed that the error rate of simulated reads is higher than QV.
    - FASTQSim:
      - Takes long time in simulating
      - Not flexibly to directly change parameters.
      - Error rate of simulated sequences produced not well matched to the QV.
  - Method:
    - Steps:
      - modeling read length distribution
      - selecting QVs
      - calculating overall base error probability
      - assigning different base error probabilities.
    - Read length:
      - log-normal distribution
      - sequencing depth
      - sampling fron known PacBio FASTA/FASTQ
    - QVs:
      - from the pre-trained QVs table
      - from user-defined table
    - Error probabilities:
      - QV of each base in read sequence is logarithmically related to the base error probability
      
### Nanosim

- Metadata:
  - Python, R
  - Oxford Nanopore MinION
  - <https://www.bcgsc.ca/resources/software/nanosim> \& <https://github.com/bcgsc/NanoSim>
  - Yang, C., Chu, J., Warren, R. L., & Birol, I. (2017). NanoSim: nanopore sequence read simulator based on statistical characterization. GigaScience, 6(4), 1–6. <https://doi.org/10.1093/gigascience/gix010>
  - Hafezqorani, S., Yang, C., Lo, T., Nip, K. M., Warren, R. L., & Birol, I. (2020). Trans-NanoSim characterizes and simulates nanopore RNA-sequencing data. GigaScience, 9(6), giaa061. <https://doi.org/10.1093/gigascience/giaa061>
  - Scopus Cited 2, 2
- Highlights:
  - pass
- Article 1:
  - A model-and-simulate approach.
  - Steps:
    - Modelling:
      - Use reference FASTA and Aligned MAF as input.
    - Simulation:
      - lengths of errors are drawn from the statistical models, and the error types are determined by a Markov chain
  - Read Length:
    - For aligned reads, typically only a middle region can be aligned, leaving the flanking head and tail regions soft-clipped from alignments.
      - The length distribution of these head and tail regions exhibits a multimodal pattern.
      - The full read length distribution can be characterized by two empirical distributions: one for the length of the aligned regions, the second for the ratio of alignment lengths to read lengths.
    - Length distributions of unaligned reads are also generated to simulate unaligned reads.
  - Sequencing Errors:
    - Stretches of substitution errors as being distributed according to Poisson distribution.
    - Indels follow Weibull distributions.
    - All error modes have a second component of geometric distribution, which we postulate describes stochastic noise.
    - substitution errors are not uniform, with a weak bias toward G and C
    - The k-mer bias of ONT reads
      - As a DNA molecule with a stretch of homopolymer sequence traverses through a nanopore, the change in electric current is not detectable or fails to be interpreted by the base-calling algorithm, leading to a deficient representation of homopolymers longer than the number of bases that can fit in the nanopores
- Article 2:
  - pass

### DeepSimulator

- Metadata:
  - Python
  - Oxford Nanopore MinION
  - <https://github.com/liyu95/DeepSimulator>
  - Li, Y., Han, R., Bi, C., Li, M., Wang, S., & Gao, X. (2018). DeepSimulator: a deep simulator for Nanopore sequencing. Bioinformatics (Oxford, England), 34(17), 2899–2908. <https://doi.org/10.1093/bioinformatics/bty223>
  - Li, Y., Wang, S., Bi, C., Qiu, Z., Li, M., & Gao, X. (2020). DeepSimulator1.5: a more powerful, quicker and lighter simulator for Nanopore sequencing. Bioinformatics (Oxford, England), 36(8), 2578–2580. <https://doi.org/10.1093/bioinformatics/btz963>
  - Scopus Cited 34, ?
- Highlights:
  - Instead of being a simulator that only mimics the result, our simulator mimics Nanopore sequencing deeply by simulating the entire processing pipeline.
  - When translating the sequences into the current signals, we build a context-dependent pore model using deep learning methods
- Article 1:
  - Review
    - SiLiCO, NanoSim, ReadSim:
      - same property of generating simulated data utilizing the input nucleotide sequence and the explicit profiles (a set of parameters, such as insertion and deletion rates, substitution rates, read lengths, error rates and quality scores)
      - Not truly capture the complex nature of the Nanopore sequencing procedure, which contains multiple stages including sample preparation, current signal collection and basecalling
  - Steps:
    - To run the simulator, the user just need to input a reference genome or assembled contigs, specifying the coverage or the number of reads.
    - The sequence would first go through a preprocessing stage, which produces several shorter sequences, satisfying the input coverage requirement and the read length distribution of real Nanopore reads.
    - Then, those sequences would pass through the signal generation module, which contains the pore model component and the signal repeating component.
      - The pore model component is used to model the expected current signal of a given k-mer (k usually equals to 5 or 6 and here we use 5-mer without loss of generality).
      - Followed by the signal repeating component to produce the simulated current signals in a mixture alpha distribution.
      - Add Gaussian noise with the user-defined variance parameter to each position of the simulated signals.
    - Finally, the simulated signal would go through [Albacore](https://community.nanoporetech.com/protocols/albacore-offline-basecalli/v/abec_2003_v1_revad_29nov2016/linux), the ONT official basecaller, to produce the final simulated reads.
  - Pore Model:
    - Currently, all the [existing pore models](https://github.com/nanoporetech/kmer_models) are context-independent, which assign each 5-mer a fixed value for the expected current signal regardless of its location on the nucleotide sequence.
    - We propose a novel context-dependent pore model, taking advantage of deep learning techniques.
  - Read Length: Not easy to distinguish, 3 distributions observed:
    - Exponential distribution to fit it (e.g. reads from the human genome).
    - Beta distribution to fit it (e.g. reads from the E. coli genome).
    - Mixture distribution with two gamma distributions to fit it (e.g. reads from the lambda phage genome).
- Article 2:
  - Review:
    - For example, though the final simulated reads have almost the same error distribution as the real reads, for some sequences, the divergence between the simulated raw signals and the real signals can be large, which can be inconvenient for the users who care about the signal outputs.
    - In addition, the Nanopore technology has evolved greatly since DS1.0 was released.
  - Updates:
    - As for the sequence generator, we updated the sample read length distribution to reflect the newest real reads’ features
    - In terms of the signal generator, we added one more pore model, the context-independent pore model, which is much faster than the previous context-dependent one.
    - We added a low-pass filter to post-process the pore model signals.
    - We added the support for the newest official basecaller, Guppy, which can support both GPU and CPU.

### nanopore_simulation

- Python
- Oxford Nanopore
- <https://github.com/crohrandt/nanopore_simulation>
- C. Rohrandt et al., "Nanopore SimulatION – a raw data simulator for Nanopore Sequencing," 2018 IEEE International Conference on Bioinformatics and Biomedicine (BIBM), 2018, pp. 1-8, doi: 10.1109/BIBM.2018.8621253.
- Scopus Cited 2

### PaSS

- Perl and C
- PacBio
- <https://cgm.sjtu.edu.cn/PaSS/>
- Zhang, W., Jia, B., & Wei, C. (2019). PaSS: a sequencing simulator for PacBio sequencing. BMC bioinformatics, 20(1), 352. <https://doi.org/10.1186/s12859-019-2901-7>
- Scopus Cited 7

### badread

- Python
- PacBio, Oxford Nanpore, named "nanopore2018", "nanopore2020", "pacbio2016"
- <https://github.com/rrwick/Badread>
- Wick, (2019). Badread: simulation of error-prone long reads. Journal of Open Source Software, 4(36), 1316, <https://doi.org/10.21105/joss.01316>
- Scopus Cited ?

### pbsim2

- C++
- PacBio CLR, Oxford Nanopore
- <https://github.com/yukiteruono/pbsim2>
- Ono, Y., Asai, K., & Hamada, M. (2021). PBSIM2: a simulator for long-read sequencers with a novel generative model of quality scores. Bioinformatics (Oxford, England), 37(5), 589–595. <https://doi.org/10.1093/bioinformatics/btaa835>
- Scopus Cited ?

### yanosim

- Metadata:
  - Python
  - Oxford Nanopore
  - <https://github.com/bartongroup/yanosim>

### loresim2

- Metadata:
  - C++
  - ?
  - <https://github.com/gt1/loresim2>
