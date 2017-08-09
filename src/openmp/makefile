CXX=g++
CXXFLAGS=-O3
INCPATH=-I.
EXECS = dflp dflp_p

all: dflp

dflp: dflp.o dflpio.o
	$(CXX) -o dflp -fopenmp dflp.o dflpio.o

dflp_dyn_rc: dflp_dyn_rc.o dflpio.o
	$(CXX) -o dflp -fopenmp dflp_dyn_rc.o dflpio.o

dflp_p: dflp_p.o dflpio_p.o
	$(CXX) -o dflp -pg -fopenmp dflp_p.o dflpio_p.o

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -fopenmp $< -o $@

%_p.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -pg -fopenmp $< -o $@

clean:
	/bin/rm -rf *.o $(EXECS) *~
