CXX=g++
CXXFLAGS=-O3
INCPATH=-I.

EXECS= dflp dflp_p


all: dflp

dflp: dflp.o dflpio.o
	$(CXX) -o dflp dflp.o dflpio.o

dflp_dyn_rc: dflp_dyn_rc.o dflpio.o
	$(CXX) -o dflp dflp_dyn_rc.o dflpio.o

dflp_p: dflp_p.o dflpio_p.o
	$(CXX) -o dflp -pg dflp_p.o dflpio_p.o

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH)  $< -o $@

%_p.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -pg  $< -o $@

clean:
	/bin/rm -rf *.o $(EXECS) *~
