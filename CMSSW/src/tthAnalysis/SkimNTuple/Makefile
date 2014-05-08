Target  = do.exe
Objects = SkimNtuple.o 

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)
CXXFLAGS = -std=c++0x

all:$(Target)

do.exe: main.cc $(Objects)
	g++ -o $@ main.cc $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

SkimNtuple.o: SkimNtuple.C SkimNtuple.h
	g++ -o $@ SkimNtuple.C -c $(ROOTFLAGS) $(CXXFLAGS)

