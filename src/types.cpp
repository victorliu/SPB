#include "SPB.hpp"
#include "util.h"
#include <cstring>

SPB::Lattice2::Lattice2(double L[4]):Lattice(2){
	Lr[0] = L[0];
	Lr[1] = L[1];
	Lr[2] = L[2];
	Lr[3] = L[3];
	// Set Lk to inverse(transpose(Lr)
	double idet = 1./(Lr[0]*Lr[3] - Lr[1]*Lr[2]);
	Lk[0] = idet * Lr[3];
	Lk[1] = idet * -Lr[1];
	Lk[2] = idet * -Lr[2];
	Lk[3] = idet * Lr[0];
}
double SPB::Lattice2::CharacteristicKLength() const{
	double len[2] = {
		hypot(Lk[0],Lk[1]),
		hypot(Lk[2],Lk[3])
	};
	return len[0] > len[1] ? len[0] : len[1];
}

SPB::Lattice3::Lattice3(double L[9]):Lattice(3){
	Lr[0] = L[0];
	Lr[1] = L[1];
	Lr[2] = L[2];
	Lr[3] = L[3];
	Lr[4] = L[4];
	Lr[5] = L[5];
	Lr[6] = L[6];
	Lr[7] = L[7];
	Lr[8] = L[8];
	// Set Lk to inverse(transpose(Lr)	
	double co[3] = {
		Lr[8]*Lr[4] - Lr[5]*Lr[7],
		Lr[5]*Lr[6] - Lr[8]*Lr[3],
		Lr[7]*Lr[3] - Lr[4]*Lr[6]
	};
	double idet = 1./(Lr[0]*co[0] + Lr[1]*co[1] + Lr[2]*co[2]);
	Lk[0] = idet * co[0];
	Lk[1] = idet * co[1];
	Lk[2] = idet * co[2];
	Lk[3] = idet * (Lr[2]*Lr[7] - Lr[8]*Lr[1]);
	Lk[4] = idet * (Lr[8]*Lr[0] - Lr[2]*Lr[6]);
	Lk[5] = idet * (Lr[1]*Lr[6] - Lr[7]*Lr[0]);
	Lk[6] = idet * (Lr[5]*Lr[1] - Lr[2]*Lr[4]);
	Lk[7] = idet * (Lr[2]*Lr[3] - Lr[5]*Lr[0]);
	Lk[8] = idet * (Lr[4]*Lr[0] - Lr[1]*Lr[3]);
}
double SPB::Lattice3::CharacteristicKLength() const{
	double len[3] = {
		hypot3(Lk[0],Lk[1],Lk[2]),
		hypot3(Lk[3],Lk[4],Lk[5]),
		hypot3(Lk[6],Lk[7],Lk[8])
	};
	return max3(len[0],len[1],len[2]);
}













SPB::LorentzPole::LorentzPole():
	omega_0(0),
	Gamma(0),
	omega_p(0)
{
}

SPB::LorentzPole::LorentzPole(const LorentzPole &lp):
	omega_0(lp.omega_0),
	Gamma(lp.Gamma),
	omega_p(lp.omega_p)
{
}








SPB::ConstitutiveTensor::ConstitutiveTensor():
	type(SPB::ConstitutiveTensor::SCALAR)
{
	value[0] = 1.;
}
SPB::ConstitutiveTensor::ConstitutiveTensor(const ConstitutiveTensor &et):
	type(et.type)
{
	memcpy(value, et.value, sizeof(SPB::complex_t) * 9);
}



SPB::Material::Material(
	const std::string &name,
	const ConstitutiveTensor &eps,
	const std::vector<LorentzPole> &poles
):
	name(name),
	eps_inf(eps),
	poles(poles)
{
}
SPB::Material::Material(const Material &mat):
	name(mat.name),
	eps_inf(mat.eps_inf),
	poles(mat.poles)
{
}