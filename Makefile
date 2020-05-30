# unique-distances/Makefile

unique-distances.exe: unique-distances.cc
	g++ -g -Wall -std=c++11 -o $@ $<

# EOF
