export CXXFLAGS ?= -w -O3 -static -static-libgcc -static-libstdc++
export DOCKER_IMPL ?= docker

all: trinity.sif yasim.sif

# Trinity for RNA-Seq assembly

# Using Trinity from Broad Institute
trinity.sif:
	wget https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/archived/trinityrnaseq.v2.15.0.simg -O trinity.sif



.PHONY: pbsim_binaries
pbsim_binaries:
	$(MAKE) -C pbsim_builder binaries

yasim-latest.d:
	rm -fr yasim-latest.d
	git clone ../ yasim-latest.d

yasim.$(DOCKER_IMPL).timestamp: yasim-latest.d
	$(DOCKER_IMPL) build . -f yasim.Dockerfile -t yasim:latest
	touch yasim.$(DOCKER_IMPL).timestamp




