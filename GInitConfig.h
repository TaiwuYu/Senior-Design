////////////////////////
// GInitConfig.h	////
// Initial Configuration
// Y.Gao 05-20-2009	////
////////////////////////

#ifndef GInitConfig_H
#define GInitConfig_H
#define Pi 3.1416

#include "iostream"
#include "fstream"
#include "math.h"
//#include "Gmt.h"
using namespace std;

enum init_option {HOMO,SINGLE,TWIN_13,USER,LOAD,GRAIN_3,GRAIN_4,DOUBLE,TWIN_12,TWIN_159,LAM_123,LOAD_3264C,TWO,HPHASE,USER_SET,RANDOM_PARTICLE};

class GInitConfig
{
public:
	GInitConfig();
	~GInitConfig();
	
//	int Set(int n[],init_option op,
//			float *eta[],float a[],char eta_file[],
//			int *gs,char gs_file[],char input_eta[]=NULL); // parameters initialization
	int Set(int n[],init_option op,
			float *eta[],float a[],char eta_file[],
			int *gs,char gs_file[],float *conc1,float *conc2); // parameters initialization
	int Set_n(int n[]); // Set {ng,nv,nx,ny,nz}

	int Init_HOMO(float *eta[]);
	int Init_SINGLE(float *eta[],float a[],int *gs);
	int Init_TWIN_13(float *eta[],float a[],int *gs);
	int Init_USER(float *eta[],float a[],int *gs,float *conc1,float *conc2);
	int Init_LOAD(float *eta[],char eta_file[]);
	int Init_GRAIN_3(float *eta[]);
	int Init_GRAIN_4(float *eta[]);
	int Init_DOUBLE(float *eta[],float a[],int *gs);
	int Init_TWIN_12(float *eta[],float a[],int *gs);
	int Init_TWIN_159(float *eta[],float a[],int *gs);
	int Init_LAM_123(float *eta[],char eta_file[],float a[],int *gs);
	int Init_LOAD_3264C(float *eta[],char eta_file[]);
	int Init_TWO(float *eta[],float a[],int *gs);
	int Init_Hphase(float *eta[],float a[],int *gs);
       // void Init_field(float *eta[],char input_eta[]);
        void Init_field(float *eta[],float *conc1,float *conc2);
        void init_from_file(float *tpr,char *filename);
	int Init_particle(float *eta[],float a[],int *gs);
	int Output_VTK_header(ofstream*,int,int,int);


	int Init_gs(int *gs,char gs_file[]);

	inline long int index1(int g1,int v1){	return g1*nv+v1;	}
	inline long int index2(int i,int j,int k){	return (i*ny+j)*nz+k;	}

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

public:
	int ng,nv,nx,ny,nz; // number of grains, variants, x, y, z
	unsigned long int random_seed;
	float	random_ave,random_var;

	bool initflag;
	init_option option; // HOMO, SINGLE, USER, LOAD
};

#endif

