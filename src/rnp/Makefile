include ../Makefile.inc

CXXFLAGS += $(RNP_INC)

libSPBrnp: \
	z-Eigensystems.o \
	z-GeneralizedEigensystems.o \
	z-GMRES.o \
	z-IRA.o \
	z-JDQZ.o
	$(AR) rcs libSPBrnp.a z-*.o

z-Eigensystems.o: Eigensystems.cpp
	$(CXX) $(CXXFLAGS) -c Eigensystems.cpp -o z-Eigensystems.o
z-GeneralizedEigensystems.o: GeneralizedEigensystems.cpp
	$(CXX) $(CXXFLAGS) -c GeneralizedEigensystems.cpp -o z-GeneralizedEigensystems.o
z-GMRES.o: GMRES.cpp GMRES.h
	$(CXX) $(CXXFLAGS) -c GMRES.cpp -o z-GMRES.o
z-IRA.o: IRA.cpp
	$(CXX) $(CXXFLAGS) -c IRA.cpp -o z-IRA.o
z-JDQZ.o: JDQZ.cpp JDQZ.h GMRES.h
	$(CXX) $(CXXFLAGS) -c JDQZ.cpp -o z-JDQZ.o

all:
	libSPBrnp

clean:
	rm -f z-*.o
