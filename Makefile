CXX := g++
CXXFLAGS := -std=c++11 -fopenmp -Ofast -Wall -Wextra -pedantic

all: REPrise REPrise_g

REPrise: REPrise_v006.cpp sais_long.c cmd_line_opts.c
	$(CXX) $(CXXFLAGS) $^ -o $@	
	
REPrise_g: REPrise_v006.cpp sais_long.c cmd_line_opts.c
	$(CXX) $(CXXFLAGS) -g $^ -o $@

clean:
	rm REPrise REPrise_g