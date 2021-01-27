////////////////////////
// Gelastic.h	////////
// Elastic Energy	////
// Y.Gao 05-22-2009	////
////////////////////////

#ifndef Gelastic_H
#define Gelastic_H

#include "iostream"
#include "assert.h"
#include "fstream"
#include <stdlib.h>

#include "Gfftw3.h"

using namespace std;

#define Pi 3.1416
#define TINY (1E-20)
#define RSQ2 0.70710678118654744
#define C6_loop for(int mi=0;mi<6;mi++)for(int mj=0;mj<6;mj++)
#define C4_loop for(int mi=0;mi<3;mi++)for(int mj=0;mj<3;mj++)for(int mk=0;mk<3;mk++)for(int ml=0;ml<3;ml++)

enum boundary_condition {FIX,RELAX,CONFIG_DEPT};
enum elastic_couple {PHI1,PHI2,PHI23,PHI345};

class Gelastic
{
public:
	Gelastic();
	~Gelastic();
        typedef double real;
      //  typedef cufftDoubleComplex complx;
        typedef real ten2nd[3][3];
      //  typedef complx ten2ndk[3][3];
        typedef real ten4th[3][3][3][3];
        typedef real voigt[6];
       // typedef complx voigtk[6];
        typedef real voigt66[6][6];

	int Init_Space();
	int Destroy_Space();
	int Destroy_Bpq();
	int Set(int n[],float cc[],float aa[],
			float **rott,float s_app_magg,float s_app_inc,float e_app_mag,float scalee,
			boundary_condition bcc,elastic_couple ecc,int cdd,fftw_option op,
			int output_s,int output_ie);
	int Set_n(int n[]);
	int Set_Cijkl(float cc[]);
	int Set_a(float aa[]);
	int Set_SFTS(); // Stress Free Transformation Strain
	int Set_SFTS_habit();
	int Set_SFTS_NiPtTi(); // Both chemical and structural field dependent
	int Set_SFTS_PLAMINATE(); // Laminate structure for P phase
	int Set_SFTS_TiAlV();  // Used for TiAl (test for Rongpei)
	int Set_SFTS_NiTiHf();
	
	int Set_Sijkl_Hphase();
	
	int Set_u();
	int Set_s0();
	int Set_rot(float**);
	int Set_g();

	float Energy();
	int Potential(float **eta,float **d_eta,int *gs);
	int Potential_k(fftw_complex **eta,fftw_complex **d_eta_k);
	int Calc_phi_k(float **eta);

	int Calc_Bpq();
	int Calc_Bpd();
	int Output_Stress(float **eta,int *gs);
	int Output_Strain(float **eta,int *gs, int out);
	int Output_InterEnergy();
	int Output_TranStrain(float **eta,int *gs);
	int Load_StressField(float **stress);
	int Reset_SFTS_MT(); // for InteractionEnergy calculation in NiTiPt
        void LU_dcmp(real **a, int n, int *indx, real *d);
        void LU_bksb(real **a, int n, int *indx, real b[]);
        void LU_inv_66(voigt66 c);
        

	int Output_VTK_header(ofstream*,int,int,int);
        int UpdateStress(float s_applied[][3]);
	int AppliedStress(float s_applied[][3]);
	int AppliedStrain(float e_applied[][3]);
	int UpdateStrain(float e_applied[][3]);

	
	inline long int index1(int g1,int v1){	return g1*nv+v1;	}
	inline long int index2(int i,int j,int k){	return (i*ny+j)*nz+k;	}
	inline long int index11(int g1,int v1,int g2,int v2){ return ((g1*nv+v1)*ng+g2)*nv+v2; }

//        static inline void chg_basis_Kelvin_4(ten4th T4, voigt66 C2);
//        static inline void chg_basis_Kelvin_3( voigt66 C2, ten4th T4);

public:
	
	int ng,nv,nx,ny,nz;
	float c[3]; // c11,c12,c44
	float Cijkl[3][3][3][3];
	float Sijkl[3][3][3][3];
	float a[4]; // lattice mismatch
	float *g[3],*g_sqr; // g[3].[nxyz]
	float **u,**e0,**s0; // TransMatrix,SFTS,Strs [g1][v1].[3][3]
	float **rot; //rotation matrix
	float **Bpq; //[g1][v1][g2][v2].[nxyz]
	float **Bpd[3][3]; //[mm][nn].[g1][v1].[nxyz]
	fftw_complex **phi_k; //[g1][v1].[nxyz]
	float s_app[3][3],*s_app_miu,s_app_mag,s_app_inc; // applied stress,potential,magnitude
	float e_app[3][3],*e_app_miu,e_app_mag,e_app_inc;
	float **s_field;
	float scale;
	bool initflag;
	boundary_condition bc;
	elastic_couple ec;
	int cd; // config-dept
	Gfftw3 gfft;
	int output_stress;
	int output_interenergy;
};

#endif

