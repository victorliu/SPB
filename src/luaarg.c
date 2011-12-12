#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "luaarg.h"

#include <stdlib.h>
#include <string.h>

#include <stdio.h>

typedef void (*luaarg_handler)(lua_State *L, int index, void *val);

static void luaarg_handler_table(lua_State *L, int index, void *val){
	luaarg_parse(L, index, (const luaarg_argspec*)val);
}
void luaarg_handler_int(lua_State *L, int index, void *val){
	*(int*)val = luaL_checkint(L, index);
}
void luaarg_handler_double(lua_State *L, int index, void *val){
	*(double*)val = luaL_checknumber(L, index);
}
void luaarg_handler_string(lua_State *L, int index, void *val){
	char **buf = (char**)val;
	*buf = strdup(luaL_checkstring(L, index));
}
void luaarg_handler_double_vec2(lua_State *L, int index, void *val){
	int i;
	double *v = (double*)val;
	luaL_argcheck(L, lua_istable(L,index) && 2==lua_objlen(L,index), index, "Expected 2-vector");
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		v[i] = luaL_checknumber(L, -1);
		lua_pop(L, 1);
	}
}
void luaarg_handler_complex(lua_State *L, int index, void *val){
	int i;
	double *v = (double*)val;
	luaL_argcheck(L, lua_isnumber(L,index) || (lua_istable(L,index) && 2==lua_objlen(L,index)), index, "Expected complex number");
	if(lua_istable(L, index)){
		for(i = 0; i < 2; ++i){
			lua_pushinteger(L, i+1);
			lua_gettable(L, -2);
			v[i] = luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
	}else{
		v[0] = luaL_checknumber(L, index);
		v[1] = 0;
	}
}
void luaarg_handler_double_vec3(lua_State *L, int index, void *val){
	int i;
	double *v = (double*)val;
	luaL_argcheck(L, lua_istable(L,index) && 3==lua_objlen(L,index), index, "Expected 3-vector");
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		v[i] = luaL_checknumber(L, -1);
		lua_pop(L, 1);
	}
}
void luaarg_handler_double_mat2(lua_State *L, int index, void *val){
	int i,j;
	double *m = (double*)val;
	luaL_argcheck(L, lua_istable(L,index) && 2==lua_objlen(L,index), index, "Expected 2x2 matrix");
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		luaL_argcheck(L, lua_istable(L,-1) && 2==lua_objlen(L,-1), index, "Expected 2x2 matrix");
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			m[i+j*2] = luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
}
void luaarg_handler_double_mat3(lua_State *L, int index, void *val){
	int i,j;
	double *m = (double*)val;
	luaL_argcheck(L, lua_istable(L,index) && 3==lua_objlen(L,index), index, "Expected 3x3 matrix");
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		luaL_argcheck(L, lua_istable(L,-1) && 3==lua_objlen(L,-1), index, "Expected 3x3 matrix");
		for(j = 0; j < 3; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			m[i+j*3] = luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
}
void luaarg_handler_int_mat2(lua_State *L, int index, void *val){
	int i,j;
	int *m = (int*)val;
	luaL_argcheck(L, lua_istable(L,index) && 2==lua_objlen(L,index), index, "Expected 2x2 integer matrix");
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		luaL_argcheck(L, lua_istable(L,-1) && 2==lua_objlen(L,-1), index, "Expected 2x2 integer matrix");
		for(j = 0; j < 2; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			m[i+j*2] = luaL_checkint(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
}
void luaarg_handler_int_mat3(lua_State *L, int index, void *val){
	int i,j;
	int *m = (int*)val;
	luaL_argcheck(L, lua_istable(L,index) && 3==lua_objlen(L,index), index, "Expected 3x3 integer matrix");
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		luaL_argcheck(L, lua_istable(L,-1) && 3==lua_objlen(L,-1), index, "Expected 3x3 integer matrix");
		for(j = 0; j < 3; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			m[i+j*3] = luaL_checkint(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
}
void luaarg_handler_int_vec2(lua_State *L, int index, void *val){
	int i;
	int *v = (int*)val;
	luaL_argcheck(L, lua_istable(L,index) && 2==lua_objlen(L,index), index, "Expected 2-vector");
	for(i = 0; i < 2; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		v[i] = luaL_checkint(L, -1);
		lua_pop(L, 1);
	}
}
void luaarg_handler_int_vec3(lua_State *L, int index, void *val){
	int i;
	int *v = (int*)val;
	luaL_argcheck(L, lua_istable(L,index) && 3==lua_objlen(L,index), index, "Expected 3-vector");
	for(i = 0; i < 3; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		v[i] = luaL_checkint(L, -1);
		lua_pop(L, 1);
	}
}
void luaarg_handler_double_vec(lua_State *L, int index, void *val){
	int i,j;
	luaarg_double_vector *v = (luaarg_double_vector*)val;
	luaL_argcheck(L, lua_istable(L,index), index, "Expected vector");
	
	v->n = lua_objlen(L,index);
	
	v->v = (double*)malloc(sizeof(double) * v->n);
	
	for(i = 0; i < v->n; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		v->v[i] = luaL_checknumber(L, -1);
		lua_pop(L, 1);
	}
}
void luaarg_handler_double_mat(lua_State *L, int index, void *val){
	int i,j;
	luaarg_double_matrix *m = (luaarg_double_matrix*)val;
	luaL_argcheck(L, lua_istable(L,index), index, "Expected matrix");
	
	m->r = lua_objlen(L,index);
	lua_pushinteger(L, 1);
	lua_gettable(L, -2);
	m->c = lua_objlen(L, -1);
	lua_pop(L, 1);
	
	m->m = (double*)malloc(sizeof(double) * m->r * m->c);
	
	for(i = 0; i < m->r; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		luaL_argcheck(L, lua_istable(L,-1) && m->c==lua_objlen(L,-1), index, "Expected matrix");
		for(j = 0; j < m->c; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			m->m[i+j*m->r] = luaL_checknumber(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
}
void luaarg_handler_int_vec(lua_State *L, int index, void *val){
	int i,j;
	luaarg_int_vector *v = (luaarg_int_vector*)val;
	luaL_argcheck(L, lua_istable(L,index), index, "Expected integer vector");
	
	v->n = lua_objlen(L,index);
	
	v->v = (int*)malloc(sizeof(int) * v->n);
	
	for(i = 0; i < v->n; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		v->v[i] = luaL_checkint(L, -1);
		lua_pop(L, 1);
	}
}
void luaarg_handler_int_mat(lua_State *L, int index, void *val){
	int i,j;
	luaarg_int_matrix *m = (luaarg_int_matrix*)val;
	luaL_argcheck(L, lua_istable(L,index), index, "Expected integer matrix");
	
	m->r = lua_objlen(L,index);
	lua_pushinteger(L, 1);
	lua_gettable(L, -2);
	m->c = lua_objlen(L, -1);
	lua_pop(L, 1);
	
	m->m = (int*)malloc(sizeof(int) * m->r * m->c);
	
	for(i = 0; i < m->r; ++i){
		lua_pushinteger(L, i+1);
		lua_gettable(L, -2);
		luaL_argcheck(L, lua_istable(L,-1) && m->c==lua_objlen(L,-1), index, "Expected matrix");
		for(j = 0; j < m->c; ++j){
			lua_pushinteger(L, j+1);
			lua_gettable(L, -2);
			m->m[i+j*m->r] = luaL_checkint(L, -1);
			lua_pop(L, 1);
		}
		lua_pop(L, 1);
	}
}

void luaarg_handler_complex_constitutive_tensor3x3(lua_State *L, int index, void *val){
	int i, j, k;
	double *v = (double*)val;
	memset(v, 0, sizeof(double)*2*3*3);
	if(lua_isnumber(L, index)){
		v[0+(0+0*3)*2] = lua_tonumber(L, index);
		v[0+(1+1*3)*2] = v[0+(0+0*3)*2];
		v[0+(2+2*3)*2] = v[0+(0+0*3)*2];
	}else if(lua_istable(L, index)){
		const int len1 = lua_objlen(L, index);
		if(2 == len1){
			lua_pushinteger(L, 1);
			lua_gettable(L, -2);
			v[0+(0+0*3)*2] = lua_tonumber(L, -1);
			lua_pop(L, 1);
			
			lua_pushinteger(L, 2);
			lua_gettable(L, -2);
			v[1+(0+0*3)*2] = lua_tonumber(L, -1);
			lua_pop(L, 1);
			
			v[0+(1+1*3)*2] = v[0+(0+0*3)*2];
			v[1+(1+1*3)*2] = v[1+(0+0*3)*2];
			v[0+(2+2*3)*2] = v[0+(0+0*3)*2];
			v[1+(2+2*3)*2] = v[1+(0+0*3)*2];
		}else if(3 == len1){
			int isreal;
			int isdiag;
			lua_pushinteger(L, 1);
			lua_gettable(L, -2);
			if(lua_isnumber(L, -1)){
				isreal = 1;
				isdiag = 1;
			}else if(lua_istable(L, -1)){
				int len2 = lua_objlen(L, -1);
				if(2 == len2){
					isreal = 0;
					isdiag = 1;
				}else if(3 == len2){
					lua_pushinteger(L, 1);
					lua_gettable(L, -2);
					if(lua_isnumber(L, -1)){
						isreal = 1;
						isdiag = 0;
					}else{
						isreal = 0;
						isdiag = 0;
					}
					lua_pop(L, 1);
				}else{
					luaL_error(L, "Expected 3x3 complex constitutive tensor\n");
				}
			}else{
				luaL_error(L, "Expected 3x3 complex constitutive tensor\n");
			}
			lua_pop(L, 1);
			
			if(isdiag){ // diagonal and real
				for(i = 0; i < 3; ++i){
					lua_pushinteger(L, i+1);
					lua_gettable(L, -2);
					if(isreal){
						v[0+(i+i*3)*2] = luaL_checknumber(L, -1);
					}else{
						if(!lua_istable(L, -1) || 2 != lua_objlen(L, -1)){
							luaL_error(L, "Expected 3x3 complex constitutive tensor\n");
						}
						for(j = 0; j < 2; ++j){
							lua_pushinteger(L, j+1);
							lua_gettable(L, -2);
							v[j+(i+i*3)*2] = luaL_checknumber(L, -1);
							lua_pop(L, 1);
						}
					}
					lua_pop(L, 1);
				}
			}else{ // tensor, or complex diagonal
				for(i = 0; i < 3; ++i){
					lua_pushinteger(L, i+1);
					lua_gettable(L, -2);
					if(!lua_istable(L, -1) || 3 != lua_objlen(L, -1)){
						luaL_error(L, "Expected 3x3 complex constitutive tensor\n");
					}
					for(j = 0; j < 3; ++j){
						lua_pushinteger(L, i+1);
						lua_gettable(L, -2);
						if(isreal){
							v[0+(i+j*3)*2] = luaL_checknumber(L, -1);
						}else{
							if(!lua_istable(L, -1) || 2 != lua_objlen(L, -1)){
								luaL_error(L, "Expected 3x3 complex constitutive tensor\n");
							}
							for(k = 0; k < 2; ++k){
								lua_pushinteger(L, k+1);
								lua_gettable(L, -2);
								v[k+(i+j*3)*2] = luaL_checknumber(L, -1);
								lua_pop(L, 1);
							}
						}
						lua_pop(L, 1);
					}
					lua_pop(L, 1);
				}
			}
		}else{
			luaL_error(L, "Expected 3x3 complex constitutive tensor\n");
		}
	}
}

void luaarg_parse(lua_State *L, int index, const luaarg_argspec *argspec){
	static const luaarg_handler handler[] = {
		&luaarg_handler_table,
		&luaarg_handler_int,
		&luaarg_handler_double,
		&luaarg_handler_string,
		&luaarg_handler_double_vec2,
		&luaarg_handler_double_vec3,
		&luaarg_handler_double_mat2,
		&luaarg_handler_double_mat3,
		&luaarg_handler_int_vec2,
		&luaarg_handler_int_vec3,
		&luaarg_handler_int_mat2,
		&luaarg_handler_int_mat3,
		&luaarg_handler_complex,
		&luaarg_handler_complex_constitutive_tensor3x3,
		&luaarg_handler_double_vec,
		&luaarg_handler_double_mat,
		&luaarg_handler_int_vec,
		&luaarg_handler_int_mat
	};
	const luaarg_argspec *arg = argspec;
	if(index < 0){ index += 1+lua_gettop(L); }
	luaL_argcheck(L, lua_istable(L,index), index, "Expected table (named arguments)");
	while(NULL != arg->name){
		lua_getfield(L, index, arg->name);
		if(lua_isnil(L, index+1) && !arg->optional){
			lua_pop(L, 1);
			luaL_error(L, "Could not find required argument: %s\n", arg->name);
		}else{
			handler[arg->type](L, index+1, arg->valptr);
			lua_pushnil(L);
			lua_setfield(L, index, arg->name);
		}
		lua_pop(L, 1);
		arg++;
	}
	
	lua_pushnil(L);
	while(lua_next(L, index) != 0){
		/* uses 'key' (at index -2) and 'value' (at index -1) */
		luaL_error(L, "Unrecognized argument: %s\n", lua_tostring(L, -2));
		/* removes 'value'; keeps 'key' for next iteration */
		lua_pop(L, 1);
	}
}
