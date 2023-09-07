.PHONY: all
all: labw_utils

labw_utils:
	git clone --branch devel https://codeup.aliyun.com/646607f7b60a10415c320130/Bio/labw_utils.git

.PHONY: update
update:
	cd labw_utils && git pull && cd ..

# axel "https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.5.tar.gz"
# tar xvzf RepeatMasker-4.1.5.tar.gz
# cp ../../explore/test_de_novo_transcriptome_assembly/ref/Dfam_curatedonly.h5 Libraries/Dfam.h5
# ./configure -trf_prgm $(which trf) -rmblast_dir $CONDA_PREFIX/bin -hmmer_dir $CONDA_PREFIX/bin
# axel https://github.com/pedronachtigall/CodAn/archive/refs/heads/master.zip
