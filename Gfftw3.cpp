/////////////////////////
//	Gfftw3.cpp	/////////
// Y.Gao 05-18-2009	/////
/////////////////////////

#include "Gfftw3.h"
#include <stdlib.h>

Gfftw3::Gfftw3()
{
	nx=ny=nz=nxyz=0;
	option=ESTIMATE;

	initflag=false;
	data=NULL;
	plan_forward=NULL;
	plan_backward=NULL;

}

Gfftw3::~Gfftw3()
{
	if(initflag==true)
		fftw_free(data);
	data=NULL;
}

int Gfftw3::Set(int n[],fftw_option op)
{
	Set_nxyz(n);
	option=op;

	data=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nxyz);
	assert(data!=NULL);

	switch(option)
	{
	case ESTIMATE:
		plan_forward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_FORWARD,FFTW_ESTIMATE); 
		plan_backward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_BACKWARD,FFTW_ESTIMATE);
		break;
	case MEASURE:
		plan_forward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_FORWARD,FFTW_MEASURE);
		plan_backward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_BACKWARD,FFTW_MEASURE);
		break;
	case PATIENT:
		plan_forward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_FORWARD,FFTW_PATIENT);
		plan_backward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_BACKWARD,FFTW_PATIENT);
		break;
	case EXHAUSTIVE:
		plan_forward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_FORWARD,FFTW_EXHAUSTIVE);
		plan_backward=fftw_plan_dft_3d(nx,ny,nz,data,data,FFTW_BACKWARD,FFTW_EXHAUSTIVE);
		break;
	default:
		cout<<"Invalid FFT option!"<<endl;
		exit(1);
	}
	initflag=true;

	return 0;
}

int Gfftw3::Set_nxyz(int n[])
{
	if(n[2]<1 || n[3]<1 || n[4]<1)
	{
		cout<<"Invalid dimensions in Gfftw3!"<<endl;
		exit(1);
	}

	nx=n[2];
	ny=n[3];
	nz=n[4];
	nxyz=nx*ny*nz;

	return 0;
}

int Gfftw3::FFTW_3D(fftw_complex *in,fftw_complex *out,fftw_direction direction)
{
	int i,j,k;
	long indx;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				data[indx][0]=in[indx][0];
				data[indx][1]=in[indx][1];
			}

	if(direction==FORWARD)
		fftw_execute(plan_forward);
	else
		fftw_execute(plan_backward);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				if(direction==FORWARD)
				{
					out[indx][0]=data[indx][0];
					out[indx][1]=data[indx][1];
				}
				else
				{			
					out[indx][0]=data[indx][0]/nxyz;
					out[indx][1]=data[indx][1]/nxyz;
				}
			}

	return 0;
}

int Gfftw3::FFTW_3D(float *in,fftw_complex *out,fftw_direction direction)
{
	int i,j,k;
	long int indx;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				data[indx][0]=in[indx];
				data[indx][1]=0;
			}

	if(direction==FORWARD)
		fftw_execute(plan_forward);
	else
	{
		cout<<"Invalid input for FFTW!"<<endl;
		exit(1);
	}

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				out[indx][0]=data[indx][0];
				out[indx][1]=data[indx][1];
			}

	return 0;
}

int Gfftw3::FFTW_3D(fftw_complex *in,float *out,fftw_direction direction)
{
	int i,j,k;
	long int indx;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				data[indx][0]=in[indx][0];
				data[indx][1]=in[indx][1];
			}

	if(direction==FORWARD)
	{
		cout<<"Invalid input for FFTW!"<<endl;
		exit(1);
	}
	else
		fftw_execute(plan_backward);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				out[indx]=data[indx][0]/nxyz;
			}

	return 0;
}

int Gfftw3::Output_data()
{
	int i,j,k;
	long int indx;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx=index(i,j,k);
				cout<<data[indx][0]<<'\t'<<data[indx][1]<<endl;
			}
	
	return 0;
}

