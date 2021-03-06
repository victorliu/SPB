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
#include "util.h"
#include "mem.h"

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
static SPB_BandSolver* Lua_SPB_BandSolver_this(lua_State *L){
	SPB_BandSolver **pS = luaL_checkudata(L, 1, Lua_SPB_BandSolver_typename);
	luaL_argcheck(L, NULL !=  pS, 1, "Internal error: Got a NULL SPB_BandSolver** pointer");
	luaL_argcheck(L, NULL != *pS, 1, "Internal error: Got a NULL SPB_BandSolver* pointer");
	return *pS;
}

static int Lua_SPB_NewBandSolver(lua_State *L){
	luaarg_double_matrix Lr;
	char *pol = NULL;
	SPB_BandSolver **pS;
	const luaarg_argspec args[] = {
		{"Lattice"        , luaarg_type_DOUBLE_MAT, 0, &Lr},
		{"Polarization"   , luaarg_type_STRING    , 1, &pol},
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, 1, args);
	
	if(2 == Lr.r && 2 == Lr.c){
		if(NULL == pol){
			luaL_error(L, "Must specify Polarization for 2D bandsolver\n");
		}
		if('H' != pol[0] && 'E' != pol[0]){
			luaL_error(L, "Polarization must be 'E' or 'H'\n");
		}
	}else if(3 == Lr.r && 3 == Lr.c){
	}else{
		luaL_error(L, "Lattice must be a 2x2 or 3x3 matrix\n");
	}
	
	pS = (SPB_BandSolver**)lua_newuserdata(L, sizeof(SPB_BandSolver*));
	luaL_getmetatable(L, Lua_SPB_BandSolver_typename);
	lua_setmetatable(L, -2);
	if(2 == Lr.r){
		*pS = SPB_BandSolver_New(2, pol[0], Lr.m);
	}else{
		*pS = SPB_BandSolver_New(2, ' ', Lr.m);
	}
	
	/* clean up allocations made by luaarg_parse */
	SPB_Free(Lr.m, "Lr.m at " __FILE__ ":" STR(__LINE__));
	SPB_Free(pol, "pol at " __FILE__ ":" STR(__LINE__));
	return 1;
}

static int Lua_SPB_BandSolver__gc(lua_State *L){
	SPB_BandSolver_Destroy(Lua_SPB_BandSolver_this(L));
	return 0;
}

static int Lua_SPB_BandSolver_SetOptions(lua_State *L){
	int i;
	int res[3];
	int verb;
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	int dim = SPB_BandSolver_GetDimension(S);
	luaarg_argspec args[] = {
		{"Resolution"     , luaarg_type_INT_VEC3  , 1, &res[0]},
		{"Verbosity"      , luaarg_type_INT       , 1, &verb},
		{NULL, 0, 0, NULL}
	};
	if(2 == dim){
		args[0].type = luaarg_type_INT_VEC2;
	}
	luaarg_parse(L, 2, args);
	
	for(i = 0; i < dim; ++i){
		if(res[i] <= 0){
			luaL_error(L, "Resolution must be positive numbers; element %d is %d\n", i+1, res[i]);
		}
	}
	
	SPB_BandSolver_SetResolution(S, res);
	SPB_BandSolver_SetVerbosity(S, verb);
	
	return 0;
}

typedef struct{
	double *v;
	int n, nalloc;
} lorentz_pole_adder_data;
static int lorentz_pole_adder(lua_State *L, int key, int val, void* data){
	lorentz_pole_adder_data *dat = (lorentz_pole_adder_data*)data;
	if(dat->n >= dat->nalloc){
		dat->nalloc += 8;
		dat->v = realloc(dat->v, sizeof(double)*3*dat->nalloc);
	}
	
	luaarg_argspec args[] = {
		{"ResonanceFrequency", luaarg_type_COMPLEX, 0, &dat->v[3*dat->n+0]},
		{"PlasmaFrequency"   , luaarg_type_DOUBLE,  0, &dat->v[3*dat->n+2]},
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, val, args);
	
	(dat->n)++;
	
	return 1;
}

static int Lua_SPB_BandSolver_AddMaterial(lua_State *L){
	char *name;
	SPB_ConstitutiveTensor eps;
	double eps_inf[18];
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	luaarg_list_handler lorentz_pole_handler;
	lorentz_pole_adder_data data;
	data.v = NULL;
	data.n = 0;
	data.nalloc = 0;
	lorentz_pole_handler.func = &lorentz_pole_adder;
	lorentz_pole_handler.data = &data;
	luaarg_argspec args[] = {
		{"Name", luaarg_type_STRING, 0, &name},
		{"EpsilonInfinity", luaarg_type_COMPLEX_CONSTITUTIVE_TENSOR3X3, 0, &eps_inf[0]},
		{"LorentzPoles", luaarg_type_LIST, 1, &lorentz_pole_handler},
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, 2, args);
	
	eps.value = &eps_inf[0];
	eps.type = SPB_ConstitutiveTensor_TENSOR;
	SPB_BandSolver_AddMaterial(S, name, &eps);
	
	if(data.n > 0){
		int i;
		for(i = 0; i < data.n; ++i){
			SPB_LorentzPole pole;
			pole.omega_0 = data.v[3*i+0] * 2*M_PI;
			pole.Gamma = 2.*data.v[3*i+1];
			pole.omega_p = data.v[3*i+2] * 2*M_PI;
			SPB_BandSolver_Material_AddLorentzPole(S, name, &pole);
		}
		SPB_Free(data.v, "data.v at " __FILE__ ":" STR(__LINE__));
	}
	
	SPB_Free(name, "name at " __FILE__ ":" STR(__LINE__));
	return 0;
}
static int Lua_SPB_BandSolver_AddMaterialLorentzPole(lua_State *L){
	char *name;
	SPB_LorentzPole pole;
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	luaarg_argspec args[] = {
		{"Material", luaarg_type_STRING, 0, &name        },
		{"Omega0"  , luaarg_type_DOUBLE, 0, &pole.omega_0},
		{"Gamma"   , luaarg_type_DOUBLE, 0, &pole.Gamma  },
		{"OmegaP"  , luaarg_type_DOUBLE, 0, &pole.omega_p},
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, 2, args);
	
	SPB_BandSolver_Material_AddLorentzPole(S, name, &pole);
	
	SPB_Free(name, "name at " __FILE__ ":" STR(__LINE__));
	return 0;
}
static int Lua_SPB_BandSolver_SetRectangle(lua_State *L){
	char *name;
	double center[2];
	double hw[2];
	double angle = 0;
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	int dim = SPB_BandSolver_GetDimension(S);
	luaarg_argspec args[] = {
		{"Material"  , luaarg_type_STRING     , 0, &name     },
		{"Center"    , luaarg_type_DOUBLE_VEC2, 0, &center[0]},
		{"Halfwidths", luaarg_type_DOUBLE_VEC2, 0, &hw[0]    },
		{"Angle"     , luaarg_type_DOUBLE     , 1, &angle    },
		{NULL, 0, 0, NULL}
	};
	
	luaL_argcheck(L, 2 == dim, 1, "BandSolver dimension must be 2");
	
	luaarg_parse(L, 2, args);
	
	SPB_BandSolver_AddRectangle(S, name, center, hw, angle);
	
	SPB_Free(name, "name at " __FILE__ ":" STR(__LINE__));
	return 0;
}
static int Lua_SPB_BandSolver_OutputEpsilon(lua_State *L){
	char *filename = NULL;
	int res[3];
	char *format = NULL;
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	int dim = SPB_BandSolver_GetDimension(S);
	luaarg_argspec args[] = {
		{"Resolution", luaarg_type_INT_VEC3, 0, &res[0]  },
		{"Filename"  , luaarg_type_STRING  , 1, &filename},
		{"Format"    , luaarg_type_STRING  , 1, &format  },
		{NULL, 0, 0, NULL}
	};
	if(2 == dim){
		args[0].type = luaarg_type_INT_VEC2;
	}
	luaarg_parse(L, 2, args);
	
	SPB_BandSolver_OutputEpsilon(S, res, filename, format);
	
	SPB_Free(format, "format at " __FILE__ ":" STR(__LINE__));
	SPB_Free(filename, "filename at " __FILE__ ":" STR(__LINE__));
	return 0;
}
static int Lua_SPB_BandSolver_SetK(lua_State *L){
	unsigned i;
	double k[3];
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	unsigned dim = SPB_BandSolver_GetDimension(S);
	luaL_argcheck(L, lua_istable(L, 2) && dim == lua_objlen(L, 2), 2, "Expected k-vector");
	for(i = 0; i < dim; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		k[i] = luaL_checknumber(L, -1);
		lua_pop(L, 1);
	}
	SPB_BandSolver_SetK(S, k);
	return 0;
}
static int Lua_SPB_BandSolver_GetApproximateFrequencies(lua_State *L){
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	double range[2], tol = 1e-2;
	luaarg_argspec args[] = {
		{"TargetFrequencyRange", luaarg_type_DOUBLE_VEC2, 0, &range[0]},
		{"Tolerance"           , luaarg_type_DOUBLE     , 1, &tol     },
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, 2, args);
	
	SPB_ApproximateFrequency *list, *cur;
	int nlist;
	SPB_BandSolver_GetApproximateFrequencies(S, 2*M_PI*range[0], 2*M_PI*range[1], tol, &nlist, &list);
	cur = list;
	
	lua_createtable(L, nlist,0);
	int idx = 0;
	while(NULL != cur){
		lua_pushinteger(L, ++idx);
		lua_createtable(L, 2, 0);
			lua_pushinteger(L, 1);
			lua_createtable(L, 2, 0);
				lua_pushinteger(L, 1);
				lua_pushnumber(L, cur->lower/(2*M_PI));
				lua_settable(L, -3);
				
				lua_pushinteger(L, 2);
				lua_pushnumber(L, cur->upper/(2*M_PI));
				lua_settable(L, -3);
			lua_settable(L, -3);
			
			lua_pushinteger(L, 2);
			lua_pushinteger(L, cur->n);
			lua_settable(L, -3);
		lua_settable(L, -3);
		
		// Free the list elements as we go
		list = cur;
		cur = cur->next;
		SPB_Free(list, "list at " __FILE__ ":" STR(__LINE__));
	}
	return 1;
	/*
	if(n <= 0){
		lua_createtable(L, 0,0);
		return 1;
	}
#ifdef SPB_USING_C99_COMPLEX
	z = (SPB_complex_ptr)malloc(sizeof(double complex) * n);
#else
	z = (double*)malloc(sizeof(double) * 2*n);
#endif	
	SPB_BandSolver_GetFrequencies(S, &n, z);
	if(0 == n){ return 0; }
	
	lua_createtable(L, n, 0);
	for(i = 0; i < n; ++i){
		lua_pushinteger(L, i+1);
		lua_createtable(L, 2, 0);
#ifdef SPB_USING_C99_COMPLEX
		lua_pushinteger(L, 1);
		lua_pushnumber(L, creal(z[i]) / (2*M_PI));
		lua_settable(L, -3);
		
		lua_pushinteger(L, 2);
		lua_pushnumber(L, cimag(z[i]) / (2*M_PI));
		lua_settable(L, -3);
#else
		lua_pushinteger(L, 1);
		lua_pushnumber(L, z[2*i+0] / (2*M_PI));
		lua_settable(L, -3);
		
		lua_pushinteger(L, 2);
		lua_pushnumber(L, z[2*i+1] / (2*M_PI));
		lua_settable(L, -3);
#endif
		
		lua_settable(L, -3);
	}
	free(z);
	return 1;
	*/
}
static int Lua_SPB_BandSolver_GetBandsNear(lua_State *L){
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	double freq, tol = 1e-7;
	int n_bands;
	luaarg_argspec args[] = {
		{"TargetFrequency", luaarg_type_DOUBLE, 0, &freq   },
		{"NumBands"       , luaarg_type_INT   , 0, &n_bands},
		{"Tolerance"      , luaarg_type_DOUBLE, 1, &tol    },
		{NULL, 0, 0, NULL}
	};
	luaarg_parse(L, 2, args);
	return 0;
}
static int Lua_SPB_BandSolver_GetPerturbedFrequencies(lua_State *L){
	SPB_BandSolver *S = Lua_SPB_BandSolver_this(L);
	return 0;
}
static void Lua_SPB_Lib_Init(lua_State *L){
	static const struct luaL_Reg Lua_SPB_lib[] = {
		{"NewBandSolver", Lua_SPB_NewBandSolver},
		{NULL, NULL}
	};
	static const struct luaL_Reg Lua_SPB_BandSolver[] = {
		{"SetOptions", Lua_SPB_BandSolver_SetOptions},
		{"AddMaterial", Lua_SPB_BandSolver_AddMaterial},
		{"AddMaterialLorentzPole", Lua_SPB_BandSolver_AddMaterialLorentzPole},
		{"SetRectangle", Lua_SPB_BandSolver_SetRectangle},
		{"OutputEpsilon", Lua_SPB_BandSolver_OutputEpsilon},
		{"SetK", Lua_SPB_BandSolver_SetK},
		{"GetApproximateFrequencies", Lua_SPB_BandSolver_GetApproximateFrequencies},
		{"GetBandsNear", Lua_SPB_BandSolver_GetBandsNear},
		{"GetPerturbedFrequencies", Lua_SPB_BandSolver_GetPerturbedFrequencies},
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
	
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
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
		fprintf(stdout, "No input file given, running in interactive mode\n"); fflush(stdout);
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
