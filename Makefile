CXX := g++
CXXFLAGS := -std=c++11 -fopenmp -Ofast -Wall -Wextra -pedantic

all: REPrise

REPrise: REPrise_v1.0.0.cpp sais_long.c cmd_line_opts.c
	$(CXX) $(CXXFLAGS) $^ -o $@	

clean:
	rm REPrise