/////////////////////////
// Gfftw3.h	/////////////
// fftw3 interface	/////
// Y.Gao 05-18-2009	/////
/////////////////////////

#ifndef Gfftw3_H
#define Gfftw3_H

#include "complex"
#include "fftw3.h"
#include "iostream"
#include "assert.h"
#include <stdlib.h>

using namespace std;

enum	fftw_direction {FORWARD,BACKWARD};
enum	fftw_option {ESTIMATE,MEASURE,PATIENT,EXHAUSTIVE};

class Gfftw3
{
public:
	Gfftw3();
	~Gfftw3();

	int Set(int n[],fftw_option); // parameters initialization
	int Set_nxyz(int n[]); // Set {ng,nv,nx,ny,nz}
	int Output_data(); // Output result (complex form)

	int FFTW_3D(fftw_complex *,fftw_complex *,fftw_direction);
	int FFTW_3D(float *,fftw_complex *,fftw_direction);
	int FFTW_3D(fftw_complex *,float *,fftw_direction);

	inline long int index(int i,int j,int k) {return (i*ny+j)*nz+k;}

public:
	int nx,ny,nz; // system size
	long int nxyz; // nxyz=nx*ny*nz
	bool initflag;
	fftw_plan plan_forward;
	fftw_plan plan_backward;
	fftw_option option;
	fftw_complex* data;
};

#endif

