////////////////////////
// Gchemical.h	////////
// Chemical free energy
// Y.Gao 05-18-2009	////
////////////////////////

#ifndef Gchemical_H
#define Gchemical_H

#include "iostream"
#include "fstream"
#include "assert.h"
#include "math.h"
#include <stdlib.h>

using namespace std;

enum chemical_option {P234,P246,PREP}; // Polynomail 2,3,4 order...

class Gchemical
{
public:
	~Gchemical();
	Gchemical();

	int Set(int n[],float *aa,chemical_option);	// parameters initialization
	int Set_n(int n[]); // Set {ng,nv,nx,ny,nz}
	int Set_a(float *aa); // Set a[2]

	float Energy(float* eta[],int *gs);
	float E_P234(float* eta[],int *gs);
	float E_P246(float* eta[],int *gs);
	float E_PREP(float* eta[],int *gs);

	int Potential(float* eta[],float* d_eta[],int *gs,float *conc1,float *conc2);
	int Miu_P234(float* eta[],float* d_eta[],int *gs);
	int Miu_P246(float* eta[],float* d_eta[],int *gs,float *conc1,float *conc2);
	int Miu_PREP(float* eta[],float* d_eta[],int *gs);
	int Output_VTK_header(ofstream*,int,int,int);
        int UpdateT(float *aa); 	
        int UpdateT2(float *aa); 	
	inline long int index1(int g1,int v1){	return g1*nv+v1;	}
	inline long int index2(int i,int j,int k){	return (i*ny+j)*nz+k;	}

public:
	int ng,nv,nx,ny,nz; // number of grains, variants, x, y, z
	float a[10]; // parameters of Landau polynomial
	bool initflag;
	chemical_option option; // P234,P246
};

#endif

