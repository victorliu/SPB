include Makefile.inc

SUBDIRS = geom rnp
OBJLIBS = geom/libSPBgeom.a rnp/libSPBrnp.a

CFLAGS += -Igeom -Irnp
CXXFLAGS += -Igeom -Irnp $(RNP_INC) $(SLU_INC)

SPB: $(SUBDIRS) \
	z-main.o \
	z-luaarg.o \
	z-NestedDissection.o \
	z-LDL2.o \
	z-SPBGlue.o \
	z-util.o \
	z-types.o \
	z-IntervalEigensolver.o \
	z-BandSolver.o \
	z-BandSolver_Ez.o
	$(CXX) $(LDFLAGS) z-*.o $(OBJLIBS) $(LIBS) -o SPB

z-main.o: main.c SPB.h luaarg.h
	$(CC) $(CFLAGS) -c $(LUA_INC) $(MPI_INC) main.c -o z-main.o
z-luaarg.o: luaarg.c luaarg.h
	$(CC) $(CFLAGS) -c $(LUA_INC) luaarg.c -o z-luaarg.o
z-NestedDissection.o: NestedDissection.c NestedDissection.h
	$(CC) $(CFLAGS) -c NestedDissection.c -o z-NestedDissection.o
z-LDL2.o: LDL2.cpp LDL2.h
	$(CXX) $(CFLAGS) -c LDL2.cpp -o z-LDL2.o
z-SPBGlue.o: SPBGlue.cpp SPB.h SPB.hpp
	$(CXX) $(CXXFLAGS) -c SPBGlue.cpp -o z-SPBGlue.o
z-types.o: types.cpp SPB.hpp
	$(CXX) $(CXXFLAGS) -c types.cpp -o z-types.o
z-util.o: util.cpp util.h
	$(CXX) $(CXXFLAGS) -c util.cpp -o z-util.o
z-IntervalEigensolver.o: IntervalEigensolver.cpp IntervalEigensolver.h
	$(CXX) $(CXXFLAGS) -c IntervalEigensolver.cpp -o z-IntervalEigensolver.o
z-BandSolver.o: BandSolver.cpp SPB.hpp
	$(CXX) $(CXXFLAGS) -c BandSolver.cpp -o z-BandSolver.o
z-BandSolver_Ez.o: BandSolver_Ez.cpp SPB.hpp
	$(CXX) $(CXXFLAGS) -c BandSolver_Ez.cpp -o z-BandSolver_Ez.o
z-EigenSolver.o: EigenSolver.cpp SPB.hpp
	$(CXX) $(CXXFLAGS) -c EigenSolver.cpp -o z-EigenSolver.o
z-EigenSolver_IRA.o: EigenSolver_IRA.cpp SPB.hpp
	$(CXX) $(CXXFLAGS) -c EigenSolver_IRA.cpp -o z-EigenSolver_IRA.o

all: SPB

.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

clean:
	rm -f z-*.o SPB
