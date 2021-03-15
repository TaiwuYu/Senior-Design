////////////////////////
// Gmt.h	////////////
// Martensitic Transition
// Y.Gao 05-20-2009	////
////////////////////////

#ifndef Gmt_H
#define Gmt_H

#define Pi 3.1416

#include "iostream"
#include "fstream"
#include "assert.h"
#include "math.h"
#include "time.h"
#include "stdlib.h"
#include "string.h"

#include "Gfftw3.h"
#include "Gchemical.h"
#include "Ggradient.h"
#include "Gelastic.h"
#include "GInitConfig.h"

using namespace std;

class Gmt
{
public:
	Gmt(char*);
	~Gmt();

	int Load_Input(char*);
	int Evolution();
	int Evolution_RK2();
	int Evolution_Spectrum();
	int Init_Space(); // Allocate memory
	int Destroy_Space(); // Free memory
	int Output_eta(int out);
	int Output_VTK_header(ofstream*,int,int,int); // output vtk file header
	int Output_gs();

	int Clear_d_eta();
	int Clear_d_eta_k();
	int Clear_eta(float **);
	int Set_g();
	bool Check_gs();

	
	int LangevinNoise(int time);
	int CutLargeEta();
	
	inline float Gauss_rand() // Gauss distribution random generator
	{
		float	dRan1=rand2()*1.0/2147483647; 
		float	dRan2=rand2()*1.0/2147483647; 
		float	x1=sqrt(-2*log(dRan1))*cos(2*Pi*dRan2)*random_var+random_ave; 
		float	x2=sqrt(-2*log(dRan1))*sin(2*Pi*dRan2)*random_var+random_ave; 
		
		return x1;
	}

	inline unsigned long int rand2() // period is 2^31, should be larger than nx*ny*nz 
	{	
		unsigned long int rem=0;
		unsigned long int tp1,tp2;
		
		tp1=random_seed/127773;
		tp2=random_seed%127773;
		random_seed=tp2*16807;
		rem=2836*tp1;

		if(random_seed>rem)
			random_seed-=rem;
		else
			random_seed+=(2147483647-rem);

		return random_seed;
	}
        double unirand();
	
	inline long int index1(int g1,int v1) { return g1*nv+v1; }
	inline long int index2(int i,int j,int k) { return (i*ny+j)*nz+k; }
	inline long int index3(int i,int j,int k) { return (i*(ny+2)+j)*(nz+2)+k; }

public:
	int ng,nv,nx,ny,nz; // number of grains/variants/x/y/z
	int n[5]; //{ng,nv,nx,ny,nz}
	
	float **eta_new,**eta_old,**d_eta; // [ng*nv][nxyz],[ng*nv][nxyz],[nv][nxyz]
	float *conc1,*conc2;
	fftw_complex **eta_new_k,**eta_old_k,**d_eta_k; //[ng*nv][nxyz],[ng*nv][nxyz],[ng*nv][nxyz]
	float *g[3],*g_sqr; // g[3].[nxyz],[nxyz]
	int *gs; // [nxyz]
        float *eta_sum;//[nxyz]

	long int total_step,output_step;
	float time_step;

	unsigned long int random_seed;
	int	random_ave,random_var;
	int cut_off_step1,cut_off_step2;
	float random_coef1,random_coef2;
	int fluc_space;

	int output_stress;
	int output_interenergy;

	Gchemical chem;
	int chem_op; // P234,P246
	float chem_a[10]; // a2,a3,a4,a6
	
	Ggradient grad;
	int grad_op; // ISO_P7,ISO_P27,ANI_P7,ANI_P27
	float *kappa; // gradient coefficient

	Gelastic elas;
	float elas_c[3]; // c11,c12,c44
	float elas_a[4]; // elastic mismatch: alfa, beta, gamma
	float **rotation; // grain rotation [ng][3*3]
	float s_app[3][3];
	float s_app_mag,s_app_inc;
	float e_app[3][3];
	float e_app_mag,e_app_inc;
	float elastic_scale;
	int bc; // boundary condition
	int ec; // elastic couple option
	int cd; // configuration dependent
	int fftw_op; // ESTIMATE, MEASURE, PATIENT, EXHAUSTIVE

	GInitConfig init; 
	int init_op; // HOMO, SINGLE, USER, LOAD
	float init_a[3]; // init parameters
	char input_eta[255]; // eta loading file
	char *input_gs; // gs loading file
        char *outdir=NULL;
};

#endif

