# unique-distances/Makefile

#OPTFLAGS :=
OPTFLAGS := -O2 -DNDEBUG

unique-distances.exe: unique-distances.cc
	g++ -g -Wall $(OPTFLAGS) -std=c++11 -o $@ $<

# EOF
