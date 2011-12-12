#include "AffineGeometry.h"
#include <math.h>

int vec2_normalize(vec2 *v){
	double a[2] = {fabs(v->v[0]),fabs(v->v[1])};
	double w = (a[0] > a[1] ? a[0] : a[1]);
	if(0 == w){
		return 0;
	}else{
		a[0] /= w; a[1] /= w;
		w *= sqrt(a[0]*a[0] + a[1]*a[1]);
		w = 1./w;
		v->v[0] *= w; v->v[1] *= w;
		return 1;
	}
}
int vec3_normalize(vec3 *v){
	double a[3] = {fabs(v->v[0]),fabs(v->v[1]),fabs(v->v[2])};
	double w = (a[0] > a[1] ? (a[0] > a[2] ? a[0] : a[2]) : (a[1] > a[2] ? a[1] : a[2]));
	if(0 == w){
		return 0;
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w;
		w *= sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		w = 1./w;
		v->v[0] *= w; v->v[1] *= w; v->v[2] *= w;
		return 1;
	}
}
double vec2_dot(const vec2 *u, const vec2 *v){
	return u->v[0]*v->v[0] + u->v[1]*v->v[1];
}
double vec3_dot(const vec3 *u, const vec3 *v){
	return u->v[0]*v->v[0] + u->v[1]*v->v[1] + u->v[2]*v->v[2];
}
double vec2_cross(const vec2 *u, const vec2 *v){
	return u->v[0]*v->v[1] - u->v[1]*v->v[2];
}
void vec3_cross(const vec3 *u, const vec3 *v, vec3 *result){
	result->v[0] = u->v[1]*v->v[2] - u->v[2]*v->v[1];
	result->v[1] = u->v[2]*v->v[0] - u->v[0]*v->v[2];
	result->v[2] = u->v[0]*v->v[1] - u->v[1]*v->v[0];
}
double vec2_length(const vec2 *v){
	double a[2] = {fabs(v->v[0]),fabs(v->v[1])};
	double w = (a[0] > a[1] ? a[0] : a[1]);
	if(0 == w){
		return 0;
	}else{
		a[0] /= w; a[1] /= w;
		w *= sqrt(a[0]*a[0] + a[1]*a[1]);
		return w;
	}
}
double vec3_length(const vec3 *v){
	double a[3] = {fabs(v->v[0]),fabs(v->v[1]),fabs(v->v[2])};
	double w = (a[0] > a[1] ? (a[0] > a[2] ? a[0] : a[2]) : (a[1] > a[2] ? a[1] : a[2]));
	if(0 == w){
		return 0;
	}else{
		a[0] /= w; a[1] /= w; a[2] /= w;
		w *= sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		return w;
	}
}















#ifdef HAVE_LUA


#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <stdlib.h>
#include <math.h>
#include <LuaCommon.h>








#define INDEX_START 0

static int vec3_new(lua_State *L){
	vec3 *v;
	unsigned int i, nargs;
	v = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	
	if(lua_istable(L, -2)){
		if(lua_objlen(L, -2) != 3){
			LP_error(L, "vec3 must be created with a table of 3 numbers");
			return 0;
		}else{
			for(i = 0; i < 3; ++i){
				lua_pushinteger(L, i+1);
				lua_gettable(L, -2);
				v->v[i] = luaL_checknumber(L, -1); lua_pop(L, 1);
			}
		}
	}else if(lua_isuserdata(L, -2)){
		lua_getfield(L, LUA_REGISTRYINDEX, vec3_typename);
		if(!lua_getmetatable(L, -2) || !lua_rawequal(L, -1, -2)){
			lua_pop(L, 2);
			LP_error(L, "vec3 must be created with another vec3");
			return 0;
		}
		lua_pop(L, 2);
		vec3 *c = (vec3*)lua_touserdata(L, -2);
		v->v[0] = c->v[0];
		v->v[1] = c->v[1];
		v->v[2] = c->v[2];
	}else if(lua_gettop(L) > 3){
		v->v[0] = luaL_checknumber(L, -4);
		v->v[1] = luaL_checknumber(L, -3);
		v->v[2] = luaL_checknumber(L, -2);
	}else{
		v->v[0] = 0;
		v->v[1] = 0;
		v->v[2] = 0;
	}
	
	return 1;
}

static int vec3_set(lua_State *L){
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	lua_Number n = luaL_checknumber(L, 3);
	if(lua_isnumber(L, 2)){
		lua_Integer i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			v->v[0] = n;
		}else if(INDEX_START+1 == i){
			v->v[1] = n;
		}else if(INDEX_START+2 == i){
			v->v[2] = n;
		}else{
			LP_error(L, "Invalid index: %d", i);
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			v->v[0] = n;
		}else if('y' == s[0] && '\0' == s[1]){
			v->v[1] = n;
		}else if('z' == s[0] && '\0' == s[1]){
			v->v[2] = n;
		}else{
			LP_error(L, "Invalid index: %s", s);
		}
	}
	return 0;
}
static int vec3_get(lua_State *L) {
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	if(lua_isnumber(L, 2)){
		int i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			lua_pushnumber(L, v->v[0]); return 1;
		}else if(INDEX_START+1 == i){
			lua_pushnumber(L, v->v[1]); return 1;
		}else if(INDEX_START+2 == i){
			lua_pushnumber(L, v->v[2]); return 1;
		}else{
			return 0;
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->v[0]); return 1;
		}else if('y' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->v[1]); return 1;
		}else if('z' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->v[2]); return 1;
		}else{
			return 0;
		}
	}else{
		return 0;
	}
}

static int L_vec3_normalize(lua_State *L) {
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3_normalize(v);
	return 0;
}

static int vec3_tostring(lua_State *L) {
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	lua_pushfstring(L, "{%f,%f,%f}", v->v[0], v->v[1], v->v[2]);
	return 1;
}
static int vec3_unm(lua_State *L) {
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	
	vec3 *nv = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	nv->v[0] = -v->v[0];
	nv->v[1] = -v->v[1];
	nv->v[2] = -v->v[2];
	
	return 1;
}
static int vec3_len(lua_State *L) {
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	
	lua_pushnumber(L, vec3_length(v));
	
	return 1;
}

static int vec3_add(lua_State *L) {
	vec3 *u = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3 *v = (vec3*)luaL_checkudata(L, 2, vec3_typename);
	
	vec3 *w = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	w->v[0] = u->v[0] + v->v[0];
	w->v[1] = u->v[1] + v->v[1];
	w->v[2] = u->v[2] + v->v[2];
	
	return 1;
}

static int vec3_sub(lua_State *L) {
	vec3 *u = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3 *v = (vec3*)luaL_checkudata(L, 2, vec3_typename);
	
	vec3 *w = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	w->v[0] = u->v[0] - v->v[0];
	w->v[1] = u->v[1] - v->v[1];
	w->v[2] = u->v[2] - v->v[2];
	
	return 1;
}

static int vec3_mul(lua_State *L){
	if(lua_type(L, 1) == LUA_TNUMBER){
		double s = luaL_checknumber(L, 1);
		vec3 *v = (vec3*)luaL_checkudata(L, 2, vec3_typename);
		vec3 *w = (vec3*)lua_newuserdata(L, sizeof(vec3));
		luaL_getmetatable(L, vec3_typename);
		lua_setmetatable(L, -2);
		w->v[0] = s * v->v[0];
		w->v[1] = s * v->v[1];
		w->v[2] = s * v->v[2];
		return 1;
	}else if(lua_type(L, 2) == LUA_TNUMBER){
		double s = luaL_checknumber(L, 2);
		vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
		vec3 *w = (vec3*)lua_newuserdata(L, sizeof(vec3));
		luaL_getmetatable(L, vec3_typename);
		lua_setmetatable(L, -2);
		w->v[0] = s * v->v[0];
		w->v[1] = s * v->v[1];
		w->v[2] = s * v->v[2];
		return 1;
	}else{
		vec3 *u = (vec3*)luaL_checkudata(L, 1, vec3_typename);
		vec3 *v = (vec3*)luaL_checkudata(L, 2, vec3_typename);
		lua_pushnumber(L,
			u->v[0]*v->v[0] +
			u->v[1]*v->v[1] +
			u->v[2]*v->v[2]
		);
		return 1;
	}
}

static int vec3_pow(lua_State *L){
	vec3 *u = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3 *v = (vec3*)luaL_checkudata(L, 2, vec3_typename);
	
	vec3 *w = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	
	vec3_cross(u,v,w);
	
	return 1;
}

static int vec3_eq(lua_State *L){
	vec3 *u = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3 *v = (vec3*)luaL_checkudata(L, 2, vec3_typename);
	
	int e = (
		(u->v[0] == v->v[0]) &&
		(u->v[1] == v->v[1]) &&
		(u->v[2] == v->v[2])
		);
	
	lua_pushboolean(L, e);
	
	return 1;
}

static int vec3_div(lua_State *L){
	double s = 1./luaL_checknumber(L, 2);
	vec3 *v = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3 *w = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	w->v[0] = s * v->v[0];
	w->v[1] = s * v->v[1];
	w->v[2] = s * v->v[2];
	return 1;
}
static int vec3_index(lua_State *L) {
	const char *key = luaL_checkstring(L, 2);
	lua_getmetatable(L, 1);
	lua_getfield(L, -1, key);
	if(!lua_isnil(L, -1)){ return 1; }
	lua_settop(L, 2);
	return vec3_get(L);
}

// constructs an orthonormal frame {a,b,c} such that a is parallel to the given vector
static int vec3_makeframe(lua_State *L){
	vec3 *u = (vec3*)luaL_checkudata(L, 1, vec3_typename);
	vec3 *a = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	vec3 *b = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	vec3 *c = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	
	a->v[0] = u->v[0];
	a->v[1] = u->v[1];
	a->v[2] = u->v[2];
	vec3_normalize(a);
	
	if(fabs(u->v[0]) > fabs(u->v[1])){
		double i = 1./hypot(u->v[0],u->v[2]);
		b->v[0] = -i * u->v[2];
		b->v[1] = 0;
		b->v[2] = i * u->v[0];
	}else{
		double i = 1./hypot(u->v[1],u->v[2]);
		b->v[0] = 0;
		b->v[1] = i * u->v[2];
		b->v[2] = -i * u->v[1];
	}
	c->v[0] = a->v[1]*b->v[2] - a->v[2]*b->v[1];
	c->v[1] = a->v[2]*b->v[0] - a->v[0]*b->v[2];
	c->v[2] = a->v[0]*b->v[1] - a->v[1]*b->v[0];
	
	return 3;
}

void vec3_init(lua_State *L){
	static const struct luaL_reg vec3_lib[] = {
		{"new", vec3_new},
		{NULL, NULL}
	};
	static const struct luaL_reg vec3_lib_members[] = {
		{"get", vec3_get},
		{"set", vec3_set},
		{"normalize", L_vec3_normalize},
		{"__index", vec3_index},
		{"__newindex", vec3_set},
		{"__tostring", vec3_tostring},
		{"__unm", vec3_unm},
		{"__add", vec3_add},
		{"__sub", vec3_sub},
		{"__mul", vec3_mul},
		{"__div", vec3_div},
		{"__pow", vec3_pow},
		{"__len", vec3_len},
		{"__eq", vec3_eq},
		{"makeframe", vec3_makeframe},
		{NULL, NULL}
	};
	static const struct luaL_reg vec3_lib_meta[] = {
		{"__call", vec3_new},
		{NULL, NULL}
	};
	
	luaL_newmetatable(L, vec3_typename);
	luaL_register(L, NULL, vec3_lib_members);
	luaL_register(L, "vec3", vec3_lib);
	
	lua_getglobal(L, "vec3");
	lua_newtable(L);
	luaL_register(L, NULL, vec3_lib_meta);
	lua_setmetatable(L, -2);
	lua_pop(L, 1);
}













static int vec2_new(lua_State *L){
	vec2 *v;
	unsigned int i, nargs;
	v = (vec2*)lua_newuserdata(L, sizeof(vec2));
	luaL_getmetatable(L, vec2_typename);
	lua_setmetatable(L, -2);
	
	if(lua_istable(L, -2)){
		if(lua_objlen(L, -2) != 3){
			LP_error(L, "vec2 must be created with a table of 2 numbers");
		}else{
			for(i = 0; i < 2; ++i){
				lua_pushinteger(L, i+1);
				lua_gettable(L, -2);
				v->v[i] = luaL_checknumber(L, -1); lua_pop(L, 1);
			}
		}
	}else if(lua_isuserdata(L, -2)){
		lua_getfield(L, LUA_REGISTRYINDEX, vec2_typename);
		if(!lua_getmetatable(L, -2) || !lua_rawequal(L, -1, -2)){
			lua_pop(L, 2);
			LP_error(L, "vec2 must be created with another vec2");
			return 0;
		}
		lua_pop(L, 2);
		vec2 *c = (vec2*)lua_touserdata(L, -2);
		v->v[0] = c->v[0];
		v->v[1] = c->v[1];
	}else if(lua_gettop(L) > 2){
		v->v[0] = luaL_checknumber(L, -3);
		v->v[1] = luaL_checknumber(L, -2);
	}else{
		v->v[0] = 0;
		v->v[1] = 0;
	}
	
	return 1;
}

static int vec2_set(lua_State *L){
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	lua_Number n = luaL_checknumber(L, 3);
	if(lua_isnumber(L, 2)){
		lua_Integer i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			v->v[0] = n;
		}else if(INDEX_START+1 == i){
			v->v[1] = n;
		}else{
			LP_error(L, "Invalid index: %d", i);
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			v->v[0] = n;
		}else if('y' == s[0] && '\0' == s[1]){
			v->v[1] = n;
		}else{
			LP_error(L, "Invalid index: %s", s);
		}
	}
	return 0;
}
static int vec2_get(lua_State *L) {
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	if(lua_isnumber(L, 2)){
		int i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			lua_pushnumber(L, v->v[0]); return 1;
		}else if(INDEX_START+1 == i){
			lua_pushnumber(L, v->v[1]); return 1;
		}else{
			return 0;
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->v[0]); return 1;
		}else if('y' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->v[1]); return 1;
		}else{
			return 0;
		}
	}else{
		return 0;
	}
}

static int L_vec2_normalize(lua_State *L) {
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2_normalize(v);
	return 0;
}

static int vec2_tostring(lua_State *L) {
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	lua_pushfstring(L, "{%f,%f}", v->v[0], v->v[1]);
	return 1;
}
static int vec2_unm(lua_State *L) {
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	
	vec2 *nv = (vec2*)lua_newuserdata(L, sizeof(vec2));
	luaL_getmetatable(L, vec2_typename);
	lua_setmetatable(L, -2);
	nv->v[0] = -v->v[0];
	nv->v[1] = -v->v[1];
	
	return 1;
}
static int vec2_len(lua_State *L) {
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	
	lua_pushnumber(L, vec2_length(v));
	
	return 1;
}

static int vec2_add(lua_State *L) {
	vec2 *u = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2 *v = (vec2*)luaL_checkudata(L, 2, vec2_typename);
	
	vec2 *w = (vec2*)lua_newuserdata(L, sizeof(vec2));
	luaL_getmetatable(L, vec2_typename);
	lua_setmetatable(L, -2);
	w->v[0] = u->v[0] + v->v[0];
	w->v[1] = u->v[1] + v->v[1];
	
	return 1;
}

static int vec2_sub(lua_State *L) {
	vec2 *u = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2 *v = (vec2*)luaL_checkudata(L, 2, vec2_typename);
	
	vec2 *w = (vec2*)lua_newuserdata(L, sizeof(vec2));
	luaL_getmetatable(L, vec2_typename);
	lua_setmetatable(L, -2);
	w->v[0] = u->v[0] - v->v[0];
	w->v[1] = u->v[1] - v->v[1];
	
	return 1;
}

static int vec2_mul(lua_State *L){
	if(lua_type(L, 1) == LUA_TNUMBER){
		double s = luaL_checknumber(L, 1);
		vec2 *v = (vec2*)luaL_checkudata(L, 2, vec2_typename);
		vec2 *w = (vec2*)lua_newuserdata(L, sizeof(vec2));
		luaL_getmetatable(L, vec2_typename);
		lua_setmetatable(L, -2);
		w->v[0] = s * v->v[0];
		w->v[1] = s * v->v[1];
		return 1;
	}else if(lua_type(L, 2) == LUA_TNUMBER){
		double s = luaL_checknumber(L, 2);
		vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
		vec2 *w = (vec2*)lua_newuserdata(L, sizeof(vec2));
		luaL_getmetatable(L, vec2_typename);
		lua_setmetatable(L, -2);
		w->v[0] = s * v->v[0];
		w->v[1] = s * v->v[1];
		return 1;
	}else{
		vec2 *u = (vec2*)luaL_checkudata(L, 1, vec2_typename);
		vec2 *v = (vec2*)luaL_checkudata(L, 2, vec2_typename);
		lua_pushnumber(L,
			u->v[0]*v->v[0] +
			u->v[1]*v->v[1]
		);
		return 1;
	}
}

static int vec2_pow(lua_State *L){
	vec2 *u = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2 *v = (vec2*)luaL_checkudata(L, 2, vec2_typename);
	
	lua_pushnumber(L,
			u->v[0]*v->v[1] -
			u->v[1]*v->v[0]
		);
	
	return 1;
}

static int vec2_eq(lua_State *L){
	vec2 *u = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2 *v = (vec2*)luaL_checkudata(L, 2, vec2_typename);
	
	int e = (
		(u->v[0] == v->v[0]) &&
		(u->v[1] == v->v[1])
		);
	
	lua_pushboolean(L, e);
	
	return 1;
}

static int vec2_div(lua_State *L){
	double s = 1./luaL_checknumber(L, 2);
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2 *w = (vec2*)lua_newuserdata(L, sizeof(vec2));
	luaL_getmetatable(L, vec2_typename);
	lua_setmetatable(L, -2);
	w->v[0] = s * v->v[0];
	w->v[1] = s * v->v[1];
	return 1;
}
static int vec2_index(lua_State *L) {
	const char *key = luaL_checkstring(L, 2);
	lua_getmetatable(L, 1);
	lua_getfield(L, -1, key);
	if(!lua_isnil(L, -1)){ return 1; }
	lua_settop(L, 2);
	return vec2_get(L);
}


static int vec2_perp(lua_State *L) {
	vec2 *v = (vec2*)luaL_checkudata(L, 1, vec2_typename);
	vec2 *w = (vec2*)lua_newuserdata(L, sizeof(vec2));
	luaL_getmetatable(L, vec2_typename);
	lua_setmetatable(L, -2);
	w->v[0] = -v->v[1];
	w->v[1] =  v->v[0];
	return 1;
}



void vec2_init(lua_State *L){
	static const struct luaL_reg vec2_lib[] = {
		{"new", vec2_new},
		{NULL, NULL}
	};
	static const struct luaL_reg vec2_lib_members[] = {
		{"get", vec2_get},
		{"set", vec2_set},
		{"normalize", L_vec2_normalize},
		{"__index", vec2_index},
		{"__newindex", vec2_set},
		{"__tostring", vec2_tostring},
		{"__unm", vec2_unm},
		{"__add", vec2_add},
		{"__sub", vec2_sub},
		{"__mul", vec2_mul},
		{"__div", vec2_div},
		{"__pow", vec2_pow},
		{"__len", vec2_len},
		{"__eq", vec2_eq},
		{"perp", vec2_perp},
		{NULL, NULL}
	};
	static const struct luaL_reg vec2_lib_meta[] = {
		{"__call", vec2_new},
		{NULL, NULL}
	};
	
	luaL_newmetatable(L, vec2_typename);
	luaL_register(L, NULL, vec2_lib_members);
	luaL_register(L, "vec2", vec2_lib);
	
	lua_getglobal(L, "vec2");
	lua_newtable(L);
	luaL_register(L, NULL, vec2_lib_meta);
	lua_setmetatable(L, -2);
	lua_pop(L, 1);
}







static int pt3_new(lua_State *L){
	pt3 *v;
	unsigned int i, nargs;
	v = (pt3*)lua_newuserdata(L, sizeof(pt3));
	luaL_getmetatable(L, pt3_typename);
	lua_setmetatable(L, -2);
	
	if(lua_istable(L, -2)){
		if(lua_objlen(L, -2) != 3){
			LP_error(L, "pt3 must be created with a table of 3 numbers");
		}else{
			for(i = 0; i < 3; ++i){
				lua_pushinteger(L, i+1);
				lua_gettable(L, -2);
				v->r[i] = luaL_checknumber(L, -1); lua_pop(L, 1);
			}
		}
	}else if(lua_isuserdata(L, -2)){
		lua_getfield(L, LUA_REGISTRYINDEX, pt3_typename);
		if(!lua_getmetatable(L, -2) || !lua_rawequal(L, -1, -2)){
			lua_pop(L, 2);
			LP_error(L, "pt3 must be created with another pt3");
			return 0;
		}
		lua_pop(L, 2);
		pt3 *c = (pt3*)lua_touserdata(L, -2);
		v->r[0] = c->r[0];
		v->r[1] = c->r[1];
		v->r[2] = c->r[2];
	}else if(lua_gettop(L) > 3){
		v->r[0] = luaL_checknumber(L, -4);
		v->r[1] = luaL_checknumber(L, -3);
		v->r[2] = luaL_checknumber(L, -2);
	}else{
		v->r[0] = 0;
		v->r[1] = 0;
		v->r[2] = 0;
	}
	
	return 1;
}

static int pt3_set(lua_State *L){
	pt3 *v = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	lua_Number n = luaL_checknumber(L, 3);
	if(lua_isnumber(L, 2)){
		lua_Integer i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			v->r[0] = n;
		}else if(INDEX_START+1 == i){
			v->r[1] = n;
		}else if(INDEX_START+2 == i){
			v->r[2] = n;
		}else{
			LP_error(L, "Invalid index: %d", i);
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			v->r[0] = n;
		}else if('y' == s[0] && '\0' == s[1]){
			v->r[1] = n;
		}else if('z' == s[0] && '\0' == s[1]){
			v->r[2] = n;
		}else{
			LP_error(L, "Invalid index: %s", s);
		}
	}
	return 0;
}
static int pt3_get(lua_State *L) {
	pt3 *v = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	if(lua_isnumber(L, 2)){
		int i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			lua_pushnumber(L, v->r[0]); return 1;
		}else if(INDEX_START+1 == i){
			lua_pushnumber(L, v->r[1]); return 1;
		}else if(INDEX_START+2 == i){
			lua_pushnumber(L, v->r[2]); return 1;
		}else{
			return 0;
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->r[0]); return 1;
		}else if('y' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->r[1]); return 1;
		}else if('z' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->r[2]); return 1;
		}else{
			return 0;
		}
	}else{
		return 0;
	}
}

static int pt3_tostring(lua_State *L) {
	pt3 *v = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	lua_pushfstring(L, "{%f,%f,%f}", v->r[0], v->r[1], v->r[2]);
	return 1;
}

static int pt3_add(lua_State *L) {
	pt3 *u = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	pt3 *v = (pt3*)luaL_checkudata(L, 2, pt3_typename);
	
	pt3 *w = (pt3*)lua_newuserdata(L, sizeof(pt3));
	luaL_getmetatable(L, pt3_typename);
	lua_setmetatable(L, -2);
	w->r[0] = u->r[0] + v->r[0];
	w->r[1] = u->r[1] + v->r[1];
	w->r[2] = u->r[2] + v->r[2];
	
	return 1;
}

static int pt3_sub(lua_State *L) {
	pt3 *u = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	pt3 *w;
	vec3 *v;
	void *p = lua_touserdata(L, 2);
	if(NULL == p || !lua_getmetatable(L, 2)){
		LP_error(L, "Expected pt3 or vec3 as second argument to subtraction");
		return 0;
	}
	
	lua_getfield(L, LUA_REGISTRYINDEX, pt3_typename);
	if(!lua_rawequal(L, -1, -2)){
		v = (vec3*)p;
		lua_pop(L, 1); // pop off pt3 typename
		lua_getfield(L, LUA_REGISTRYINDEX, vec3_typename);
		if(!lua_rawequal(L, -1, -2)){
			lua_pop(L, 2);
			LP_error(L, "Expected pt3 or vec3 as second argument to subtraction");
			return 0;
		}
		lua_pop(L, 2);
		// pt3 - vec3
		w = (pt3*)lua_newuserdata(L, sizeof(pt3));
		luaL_getmetatable(L, pt3_typename);
		lua_setmetatable(L, -2);
		w->r[0] = u->r[0] - v->v[0];
		w->r[1] = u->r[1] - v->v[1];
		w->r[2] = u->r[2] - v->v[2];
		
		return 1;
	}
	lua_pop(L, 2);
	// pt3 - pt3
	w = (pt3*)p;
	v = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	v->v[0] = u->r[0] - w->r[0];
	v->v[1] = u->r[1] - w->r[1];
	v->v[2] = u->r[2] - w->r[2];
	
	return 1;
}

static int pt3_eq(lua_State *L){
	pt3 *u = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	pt3 *v = (pt3*)luaL_checkudata(L, 2, pt3_typename);
	
	int e = (
		(u->r[0] == v->r[0]) &&
		(u->r[1] == v->r[1]) &&
		(u->r[2] == v->r[2])
		);
	
	lua_pushboolean(L, e);
	
	return 1;
}

static int pt3_index(lua_State *L) {
	const char *key = luaL_checkstring(L, 2);
	lua_getmetatable(L, 1);
	lua_getfield(L, -1, key);
	if(!lua_isnil(L, -1)){ return 1; }
	lua_settop(L, 2);
	return pt3_get(L);
}

static int pt3_interp(lua_State *L) {
	pt3 *u = (pt3*)luaL_checkudata(L, 1, pt3_typename);
	pt3 *v = (pt3*)luaL_checkudata(L, 2, pt3_typename);
	double t = luaL_checknumber(L, 3);
	
	pt3 *w = (pt3*)lua_newuserdata(L, sizeof(pt3));
	luaL_getmetatable(L, pt3_typename);
	lua_setmetatable(L, -2);
	
	w->r[0] = t*v->r[0] + (1.-t)*u->r[0];
	w->r[1] = t*v->r[1] + (1.-t)*u->r[1];
	w->r[2] = t*v->r[2] + (1.-t)*u->r[2];
	return 1;
}

static int pt3_origin(lua_State *L) {
	pt3 *w = (pt3*)lua_newuserdata(L, sizeof(pt3));
	luaL_getmetatable(L, pt3_typename);
	lua_setmetatable(L, -2);
	
	w->r[0] = 0;
	w->r[1] = 0;
	w->r[2] = 0;
	return 1;
}

void pt3_init(lua_State *L){
	static const struct luaL_reg pt3_lib[] = {
		{"new", pt3_new},
		{"origin", pt3_origin},
		{"interp", pt3_interp},
		{NULL, NULL}
	};
	static const struct luaL_reg pt3_lib_members[] = {
		{"get", pt3_get},
		{"set", pt3_set},
		{"__index", pt3_index},
		{"__newindex", pt3_set},
		{"__tostring", pt3_tostring},
		{"__add", pt3_add}, // pt + vec
		{"__sub", pt3_sub}, // pt - pt
		{"__eq", pt3_eq},
		{NULL, NULL}
	};
	static const struct luaL_reg pt3_lib_meta[] = {
		{"__call", pt3_new},
		{NULL, NULL}
	};
	
	luaL_newmetatable(L, pt3_typename);
	luaL_register(L, NULL, pt3_lib_members);
	luaL_register(L, "pt3", pt3_lib);
	
	lua_getglobal(L, "pt3");
	lua_newtable(L);
	luaL_register(L, NULL, pt3_lib_meta);
	lua_setmetatable(L, -2);
	lua_pop(L, 1);
}







static int pt2_new(lua_State *L){
	pt2 *v;
	unsigned int i, nargs;
	v = (pt2*)lua_newuserdata(L, sizeof(pt2));
	luaL_getmetatable(L, pt2_typename);
	lua_setmetatable(L, -2);
	
	if(lua_istable(L, -2)){
		if(lua_objlen(L, -2) != 2){
			LP_error(L, "pt2 must be created with a table of 3 numbers");
		}else{
			for(i = 0; i < 2; ++i){
				lua_pushinteger(L, i+1);
				lua_gettable(L, -2);
				v->r[i] = luaL_checknumber(L, -1); lua_pop(L, 1);
			}
		}
	}else if(lua_isuserdata(L, -2)){
		lua_getfield(L, LUA_REGISTRYINDEX, pt2_typename);
		if(!lua_getmetatable(L, -2) || !lua_rawequal(L, -1, -2)){
			lua_pop(L, 2);
			LP_error(L, "pt2 must be created with another pt2");
			return 0;
		}
		lua_pop(L, 2);
		pt2 *c = (pt2*)lua_touserdata(L, -2);
		v->r[0] = c->r[0];
		v->r[1] = c->r[1];
	}else if(lua_gettop(L) > 3){
		v->r[0] = luaL_checknumber(L, -3);
		v->r[1] = luaL_checknumber(L, -2);
	}else{
		v->r[0] = 0;
		v->r[1] = 0;
	}
	
	return 1;
}

static int pt2_set(lua_State *L){
	pt2 *v = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	lua_Number n = luaL_checknumber(L, 3);
	if(lua_isnumber(L, 2)){
		lua_Integer i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			v->r[0] = n;
		}else if(INDEX_START+1 == i){
			v->r[1] = n;
		}else{
			LP_error(L, "Invalid index: %d", i);
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			v->r[0] = n;
		}else if('y' == s[0] && '\0' == s[1]){
			v->r[1] = n;
		}else{
			LP_error(L, "Invalid index: %s", s);
		}
	}
	return 0;
}
static int pt2_get(lua_State *L) {
	pt2 *v = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	if(lua_isnumber(L, 2)){
		int i = lua_tointeger(L, 2);
		if(INDEX_START+0 == i){
			lua_pushnumber(L, v->r[0]); return 1;
		}else if(INDEX_START+1 == i){
			lua_pushnumber(L, v->r[1]); return 1;
		}else{
			return 0;
		}
	}else if(lua_isstring(L, 2)){
		const char *s = lua_tostring(L, 2);
		if('x' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->r[0]); return 1;
		}else if('y' == s[0] && '\0' == s[1]){
			lua_pushnumber(L, v->r[1]); return 1;
		}else{
			return 0;
		}
	}else{
		return 0;
	}
}

static int pt2_tostring(lua_State *L) {
	pt2 *v = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	lua_pushfstring(L, "{%f,%f}", v->r[0], v->r[1]);
	return 1;
}

static int pt2_add(lua_State *L) {
	pt2 *u = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	pt2 *v = (pt2*)luaL_checkudata(L, 2, pt2_typename);
	
	pt2 *w = (pt2*)lua_newuserdata(L, sizeof(pt2));
	luaL_getmetatable(L, pt2_typename);
	lua_setmetatable(L, -2);
	w->r[0] = u->r[0] + v->r[0];
	w->r[1] = u->r[1] + v->r[1];
	
	return 1;
}

static int pt2_sub(lua_State *L) {
	pt2 *u = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	pt2 *w;
	vec3 *v;
	void *p = lua_touserdata(L, 2);
	if(NULL == p || !lua_getmetatable(L, 2)){
		LP_error(L, "Expected pt2 or vec3 as second argument to subtraction");
		return 0;
	}
	
	lua_getfield(L, LUA_REGISTRYINDEX, pt2_typename);
	if(!lua_rawequal(L, -1, -2)){
		v = (vec3*)p;
		lua_pop(L, 1); // pop off pt2 typename
		lua_getfield(L, LUA_REGISTRYINDEX, vec3_typename);
		if(!lua_rawequal(L, -1, -2)){
			lua_pop(L, 2);
			LP_error(L, "Expected pt2 or vec3 as second argument to subtraction");
			return 0;
		}
		lua_pop(L, 2);
		// pt2 - vec3
		w = (pt2*)lua_newuserdata(L, sizeof(pt2));
		luaL_getmetatable(L, pt2_typename);
		lua_setmetatable(L, -2);
		w->r[0] = u->r[0] - v->v[0];
		w->r[1] = u->r[1] - v->v[1];
		
		return 1;
	}
	lua_pop(L, 2);
	// pt2 - pt2
	w = (pt2*)p;
	v = (vec3*)lua_newuserdata(L, sizeof(vec3));
	luaL_getmetatable(L, vec3_typename);
	lua_setmetatable(L, -2);
	v->v[0] = u->r[0] - w->r[0];
	v->v[1] = u->r[1] - w->r[1];
	
	return 1;
}

static int pt2_eq(lua_State *L){
	pt2 *u = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	pt2 *v = (pt2*)luaL_checkudata(L, 2, pt2_typename);
	
	int e = (
		(u->r[0] == v->r[0]) &&
		(u->r[1] == v->r[1])
		);
	
	lua_pushboolean(L, e);
	
	return 1;
}

static int pt2_index(lua_State *L) {
	const char *key = luaL_checkstring(L, 2);
	lua_getmetatable(L, 1);
	lua_getfield(L, -1, key);
	if(!lua_isnil(L, -1)){ return 1; }
	lua_settop(L, 2);
	return pt2_get(L);
}

static int pt2_interp(lua_State *L) {
	pt2 *u = (pt2*)luaL_checkudata(L, 1, pt2_typename);
	pt2 *v = (pt2*)luaL_checkudata(L, 2, pt2_typename);
	double t = luaL_checknumber(L, 3);
	
	pt2 *w = (pt2*)lua_newuserdata(L, sizeof(pt2));
	luaL_getmetatable(L, pt2_typename);
	lua_setmetatable(L, -2);
	
	w->r[0] = t*v->r[0] + (1.-t)*u->r[0];
	w->r[1] = t*v->r[1] + (1.-t)*u->r[1];
	return 1;
}

static int pt2_origin(lua_State *L) {
	pt2 *w = (pt2*)lua_newuserdata(L, sizeof(pt2));
	luaL_getmetatable(L, pt2_typename);
	lua_setmetatable(L, -2);
	
	w->r[0] = 0;
	w->r[1] = 0;
	return 1;
}

void pt2_init(lua_State *L){
	static const struct luaL_reg pt2_lib[] = {
		{"new", pt2_new},
		{"origin", pt2_origin},
		{"interp", pt2_interp},
		{NULL, NULL}
	};
	static const struct luaL_reg pt2_lib_members[] = {
		{"get", pt2_get},
		{"set", pt2_set},
		{"__index", pt2_index},
		{"__newindex", pt2_set},
		{"__tostring", pt2_tostring},
		{"__add", pt2_add}, // pt + vec
		{"__sub", pt2_sub}, // pt - pt
		{"__eq", pt2_eq},
		{NULL, NULL}
	};
	static const struct luaL_reg pt2_lib_meta[] = {
		{"__call", pt2_new},
		{NULL, NULL}
	};
	
	luaL_newmetatable(L, pt2_typename);
	luaL_register(L, NULL, pt2_lib_members);
	luaL_register(L, "pt2", pt2_lib);
	
	lua_getglobal(L, "pt2");
	lua_newtable(L);
	luaL_register(L, NULL, pt2_lib_meta);
	lua_setmetatable(L, -2);
	lua_pop(L, 1);
}

void AffineGeometry_init(lua_State *L){
	vec3_init(L);
	vec2_init(L);
	pt3_init(L);
	pt2_init(L);
}
#endif // HAVE_LUA
