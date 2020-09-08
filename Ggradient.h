//////////////////////////
// Ggradient.h ////////
// Gradient Energy	//
// Y.Gao 05-18-2009	//////
//////////////////////////

#ifndef Ggradient_H
#define Ggradient_H

#include "iostream"
#include "assert.h"
#include <stdlib.h>

using namespace std;

enum gradient_option {ISO_P7,ISO_P27,ANI_P7,ANI_P27};

class Ggradient
{
public:
	Ggradient();
	~Ggradient();

	int Set(int n[],float *ka,gradient_option); // parameters initialization
	int Set_n(int n[]); // Set {ng,nv,nx,ny,nz}
	int Set_kappa(float *ka);

	float Energy(float* eta[],int *gs);
	float E_ISO_P7(float* eta[],int *gs);
	float E_ISO_P27(float* eta[],int *gs);
	float E_ANI_P7(float* eta[],int *gs);
	float E_ANI_P27(float* eta[],int *gs);

	int Potential(float* eta[],float* d_eta[],int *gs);
	int Miu_ISO_P7(float* eta[],float* d_eta[],int *gs);
	int Miu_ISO_P27(float* eta[],float* d_eta[],int *gs);
	int Miu_ANI_P7(float* eta[],float* d_eta[],int *gs);
	int Miu_ANI_P27(float* eta[],float* d_eta[],int *gs);

	inline long int index1(int g1,int v1) { return g1*nv+v1; }
	inline long int index2(int i,int j,int k) { return (i*ny+j)*nz+k; }
	inline long int index11(int g1,int v1,int g2,int v2){ return ((g1*nv+v1)*ng+g2)*nv+v2; }

public:
	int ng,nv,nx,ny,nz; // number of grains, variants, x, y, z
	float *kappa; // gradient coefficients
	bool initflag;
	gradient_option option; // ISO_P7, ISO_P27, ANI_P7, ANI_P27
};

#endif  

