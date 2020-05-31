# unique-distances/Makefile

OPTFLAGS := -O2 -DNDEBUG

all: unique-distances.exe
all: unique-distances-opt.exe

unique-distances.exe: unique-distances.cc
	g++ -g -Wall -std=c++11 -o $@ $<

unique-distances-opt.exe: unique-distances.cc
	g++ -g -Wall $(OPTFLAGS) -std=c++11 -o $@ $<

# EOF
