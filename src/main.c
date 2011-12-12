//#include "config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "luaarg.h"
#include "SPB.h"

static void usage(){
	printf("SPB [-a arg] [-h] [-t thread-count] [-v] [input-file]\n");
}
static void version(){
	printf("Stanford Photonic Bands\n");
	//printf("Version %s\n", PACKAGE_VERSION);
#ifdef HAVE_MPI
	printf("  With MPI support\n");
#endif
}

// from unistd.h
extern char   *optarg;
extern int    optind, opterr, optopt;
int getopt(int argc, char * const argv[], const char *optstring);
#include <ctype.h>

static const char interactive_key = 'k';
static void SPB_set_interactive(lua_State *L, int i){
	lua_pushlightuserdata(L, (void *)&interactive_key);
	lua_pushnumber(L, i);
	lua_settable(L, LUA_REGISTRYINDEX);
}
static int SPB_get_interactive(lua_State *L){
	int n;
	lua_pushlightuserdata(L, (void *)&interactive_key);
	lua_gettable(L, LUA_REGISTRYINDEX);
	n = lua_tointeger(L, -1);
	lua_pop(L, 1);
	return n;
}

static const char* Lua_SPB_BandSolver_typename = "SPB.BandSolver";

static int Lua_SPB_NewBandSolver(lua_State *L){
	int i;
	luaarg_double_matrix Lr;
	char *pol = NULL;
	int numbands = 0;
	luaarg_int_vector res;
	double targ[2];
	double tol;
	const luaarg_argspec args[] = {
		{"Lattice"        , luaarg_type_DOUBLE_MAT, 0, &Lr},
		{"Polarization"   , luaarg_type_STRING    , 1, &pol},
		{"NumBands"       , luaarg_type_INT       , 1, &numbands},
		{"Resolution"     , luaarg_type_INT_VEC   , 1, &res},
		{"TargetFrequency", luaarg_type_COMPLEX   , 1, &targ[0]},
		{"Tolerance"      , luaarg_type_DOUBLE    , 1, &tol},
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, 2, args);
	
	if(2 == Lr.r && 2 == Lr.c){
		if(NULL == pol){
			luaL_error(L, "Must specify Polarization for 2D bandsolver\n");
		}
		if(2 != res.n){
			luaL_error(L, "Must specify 2-vector for resolution\n");
		}
		if('H' != pol[0] && 'E' != pol[0]){
			luaL_error(L, "Polarization must be 'E' or 'H'\n");
		}
	}else if(3 == Lr.r && 3 == Lr.c){
		if(3 != res.n){
			luaL_error(L, "Must specify 3-vector for resolution\n");
		}
	}else{
		luaL_error(L, "Lattice must be a 2x2 or 3x3 matrix\n");
	}
	if(numbands < 0){
		luaL_error(L, "NumBands must >= 0\n");
	}
	for(i = 0; i < res.n; ++i){
		if(res.v[i] <= 0){
			luaL_error(L, "Resolution must be positive numbers\n");
		}
	}
	if(tol <= 0 || tol >= 1){
		luaL_error(L, "Tolerance must be in (0,1)\n");
	}
	
	free(pol);
	return 0;
}

static int Lua_SPB_BandSolver__gc(lua_State *L){
	SPB_BandSolver *S = (SPB_BandSolver*)luaL_checkudata(L, 1, Lua_SPB_BandSolver_typename);
	SPB_BandSolver_Destroy(S);
	return 0;
}

static int Lua_SPB_BandSolver_AddMaterial(lua_State *L){
}
static void Lua_SPB_Lib_Init(lua_State *L){
	static const struct luaL_Reg Lua_SPB_lib[] = {
		{"NewBandSolver", Lua_SPB_NewBandSolver},
		{NULL, NULL}
	};
	static const struct luaL_Reg Lua_SPB_BandSolver[] = {
		{"AddMaterial", Lua_SPB_BandSolver_AddMaterial},
		{NULL, NULL}
	};
	
	luaL_register(L, "SPB", Lua_SPB_lib);
	lua_pop(L, 1);
	
	luaL_newmetatable(L, Lua_SPB_BandSolver_typename);
	lua_pushvalue(L, -1);
	lua_setfield(L, -2, "__index");
	luaL_register(L, NULL, Lua_SPB_BandSolver);
	lua_pushstring(L, "__gc");
	lua_pushcfunction(L, Lua_SPB_BandSolver__gc);
	lua_settable(L, -3);
	lua_pop(L, 1);
}

int main(int argc, char *argv[]){
	char buff[256];
	int c;
	int index, error;
	char *arg = NULL;
	int mpi_size = 1, mpi_rank = 0;
	
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
	
	opterr = 0;
	while((c = getopt(argc, argv, "a:ht:v")) != -1){
		switch(c){
		case 'a':
			arg = strdup(optarg);
			break;
		case 'h':
			usage();
			return EXIT_SUCCESS;
		case 'v':
			version();
			return EXIT_SUCCESS;
		case '?':
			if('t' == optopt){
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			}else if(isprint(optopt)){
				fprintf(stderr, "Unknown option -%c.\n", optopt);
			}else{
				fprintf(stderr, "Unknown option character '\\x%x'.\n", optopt);
			}
			usage();
			return EXIT_FAILURE;
		default:
			abort();
		}
	}
	
	lua_State *L = luaL_newstate(); /* opens Lua */
	Lua_SPB_Lib_Init(L);
	luaL_openlibs(L); /* opens the standard libraries */

	// Set the argument if there is one
	if(NULL != arg){
		lua_getglobal(L, "SPB");
		lua_pushstring(L, "arg");
		lua_pushstring(L, arg);
		lua_settable(L, -3);
		lua_pop(L, 1);
	}
	
	/*
	// Would rather have this in MPI library.
	lua_getglobal(L, "SPB");
	lua_pushstring(L, "MPIRank");
	lua_pushinteger(L, mpi_rank);
	lua_settable(L, -3);
	lua_pushstring(L, "MPISize");
	lua_pushinteger(L, mpi_size);
	lua_settable(L, -3);
	lua_pop(L, 1);
	*/
	
	if(optind < argc){ // has at least 1 argument
		SPB_set_interactive(L, 0);
		
		for(index = optind; index < argc; ++index){
			error = luaL_dofile(L, argv[index]);
			if(error){
				fprintf(stderr, "%s\n", lua_tostring(L, -1));
				lua_pop(L, 1); /* pop error message from the stack */
			}
		}
	}else{ // run in REPL mode
		fprintf(stdout, "No input file given, running in interactive mode\n");
		SPB_set_interactive(L, 1);
		while(fgets(buff, sizeof(buff), stdin) != NULL){
			error = luaL_loadbuffer(L, buff, strlen(buff), "line") || lua_pcall(L, 0, 0, 0);
			if(error){
				fprintf(stderr, "%s", lua_tostring(L, -1));
				lua_pop(L, 1); /* pop error message from the stack */
			}
		}
	}
	
	lua_close(L);
	
	if(NULL != arg){ free(arg); }
#ifdef HAVE_MPI
	MPI_Finalize();
#endif
	return EXIT_SUCCESS;
}
