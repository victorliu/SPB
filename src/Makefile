include Makefile.inc

SPB: \
	z-main.o \
	z-luaarg.o \
	z-SPBGlue.o
	$(CXX) $(LDFLAGS) z-*.o $(LIBS) -o SPB

z-main.o: main.c SPB.h luaarg.h
	$(CC) $(CFLAGS) -c $(LUA_INC) $(MPI_INC) main.c -o z-main.o
z-luaarg.o: luaarg.c luaarg.h
	$(CC) $(CFLAGS) -c $(LUA_INC) luaarg.c -o z-luaarg.o
z-SPBGlue.o: SPBGlue.cpp SPB.h SPB_private.h
	$(CXX) $(CXXFLAGS) -c $(LA_INC) SPBGlue.cpp -o z-SPBGlue.o

all:
	SPB

clean:
	rm -f z-*.o SPB
