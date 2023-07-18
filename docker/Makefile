
pbsim3.cpp:
	wget https://github.com/yukiteruono/pbsim3/raw/master/src/pbsim.cpp -O pbsim3.cpp

pbsim2.cpp:
	wget https://github.com/yukiteruono/pbsim2/raw/master/src/pbsim.cpp -O pbsim2.cpp

pbsim_builder: pbsim_builder.Dockerfile
	docker build . -f pbsim_builder.Dockerfile -t pbsim_builder:latest
	touch pbsim_builder

%: %.cpp pbsim_builder
	docker run --rm -v $(PWD):$(PWD):rw pbsim_builder:latest \
		g++ -w -O3 -static -static-libgcc -static-libstdc++ $(PWD)/$< -o $(PWD)/$@

yasim: yasim.Dockerfile pbsim2 pbsim3
	docker build . -f yasim.Dockerfile -t yasim:latest
