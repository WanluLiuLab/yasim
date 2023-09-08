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

# Must be the latest development version while building.
.PHONY: labw_utils-latest.d
labw_utils-latest.d:
	$(MAKE) -C ../deps labw_utils
	rm -fr labw_utils-latest.d
	git clone --depth 0 ../deps/labw_utils labw_utils-latest.d
	rm -rf labw_utils-latest.d/.git

# Must be the latest development version while building.
.PHONY: yasim-latest.d
yasim-latest.d:
	rm -fr yasim-latest.d
	git clone --depth 0 ../ yasim-latest.d
	rm -rf yasim-latest.d/.git

yasim.$(DOCKER_IMPL).timestamp: yasim-latest.d labw_utils-latest.d
	$(DOCKER_IMPL) build . -f yasim.Dockerfile -t yasim:latest
	touch yasim.$(DOCKER_IMPL).timestamp
