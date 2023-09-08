export CXXFLAGS ?= -w -O3 -static -static-libgcc -static-libstdc++
export DOCKER_IMPL ?= docker
export PBSIM_BUILDER_DOCKER_TAG ?= pbsim_builder:latest

.PHONY: all
all: binaries sources portable_builder

# Source code

.PHONY: sources
sources: pbsim.d pbsim3.cpp pbsim2.cpp

pbsim.d:
	git clone https://salsa.debian.org/med-team/pbsim -o pbsim.d

pbsim3.cpp:
	wget https://github.com/yukiteruono/pbsim3/raw/master/src/pbsim.cpp -O pbsim3.cpp

pbsim2.cpp:
	wget https://github.com/yukiteruono/pbsim2/raw/master/src/pbsim.cpp -O pbsim2.cpp

# The CXX builder

.PHONY: portable_builder
portable_builder: pbsim_builder.sif pbsim_builder.$(DOCKER_IMPL).tar.xz

pbsim_builder.$(DOCKER_IMPL).timestamp: pbsim_builder.Dockerfile
	$(DOCKER_IMPL) build . -f pbsim_builder.Dockerfile -t $(PBSIM_BUILDER_DOCKER_TAG)

pbsim_builder.$(DOCKER_IMPL).tar.xz: pbsim_builder.$(DOCKER_IMPL).timestamp
	$(DOCKER_IMPL) save $(PBSIM_BUILDER_DOCKER_TAG) | xz -9 -T0 -vv > pbsim_builder.$(DOCKER_IMPL).tar.xz

pbsim_builder.sif: pbsim_builder.$(DOCKER_IMPL).tar.xz
	singularity build --fakeroot pbsim_builder.sif docker-archive://pbsim_builder.$(DOCKER_IMPL).tar.xz

# Binaries

.PHONY: binaries
binaries: pbsim pbsim2 pbsim3

%: %.cpp pbsim_builder.$(DOCKER_IMPL).timestamp
	$(DOCKER_IMPL) run --rm -v $(CURDIR):$(CURDIR):rw $(PBSIM_BUILDER_DOCKER_TAG) \
		g++ $(CXXFLAGS) $(CURDIR)/$< -o $(CURDIR)/$@

pbsim: pbsim_builder.$(DOCKER_IMPL).timestamp
	$(DOCKER_IMPL) run --rm -v $(CURDIR):$(CURDIR):rw $(PBSIM_BUILDER_DOCKER_TAG) \
		g++ $(CXXFLAGS) $(CURDIR)/pbsim.d/src/*.cpp -o $(CURDIR)/$@
