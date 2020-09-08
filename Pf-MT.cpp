#include "Gfftw3.h"
#include "Gchemical.h"
#include "Ggradient.h"
#include "GInitConfig.h"
#include "Gelastic.h"
#include "Gmt.h"

#include "stdlib.h"
#include "time.h"
#include "math.h"
#include "ctime"

using namespace std;

/*
int TryGfftw3()
{
	int n[5]={1,1,3,2,1};
	float a1[30]={1,2,3,4,7,8};
	fftw_complex a2[30]={0},a3[30]={0};
	Gfftw3 aaa;
	aaa.Set(n,ESTIMATE);
	aaa.FFTW_3D(a1,a2,FORWARD);
//	aaa.FFTW_3D(a2,a3,BACKWARD);
	aaa.Output_data();

	Gfftw3 bbb;
	bbb.Set(n,MEASURE);
	bbb.FFTW_3D(a1,a2,FORWARD);
	bbb.Output_data();
	
	return 0;
}

int TryGchemical()
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	
	float a[7]={(float)0.2,0,0,0,0,0,0};
	const int ng=1,nv=3;
	const int nx=2,ny=2,nz=2;
	int n[5]={ng,nv,nx,ny,nz};
	float *eta[ng*nv],*d_eta[nv];
	int gs[nx*ny*nz]={0,0,0,0,0,0,0,0};

	for(v1=0;v1<nv;v1++)
	{
		d_eta[v1]=new float[nx*ny*nz];
		for(g1=0;g1<ng;g1++)
		{		
			indx1=g1*nv+v1;
			eta[indx1]=new float[nx*ny*nz];
		}
	}

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=i*ny*nz+j*nz+k;
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=g1*nv+v1;
						(eta[indx1])[indx2]=0;
						(d_eta[v1])[indx2]=0;
						if(i<1 && j<1)
							(eta[gs[indx2]*nv+1])[indx2]=0;
					}
			}
	(eta[0])[0]=0.5;
	Gchemical aaa;
	aaa.Set(n,a,P234);
	cout<<aaa.Energy(eta,gs)<<endl;
	aaa.Potential(eta,d_eta,gs);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=i*ny*nz+j*nz+k;
				for(v1=0;v1<nv;v1++)
				{
					for(g1=0;g1<ng;g1++)
					{
						indx1=g1*nv+v1;
						cout<<(eta[indx1])[indx2]<<'\t';
					}
				}
			}
	cout<<endl<<endl;
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=i*ny*nz+j*nz+k;
				for(v1=0;v1<nv;v1++)
					cout<<(d_eta[v1])[indx2]<<'\t';
			}
					
	for(v1=0;v1<nv;v1++)
	{
		delete[] d_eta[v1];
		for(g1=0;g1<ng;g1++)
		{		
			indx1=g1*nv+v1;
			delete[] eta[indx1];
		}
	}

	return 0;
}

int TryGgradient()
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	
	float ka[1]={1};
	const int ng=1,nv=1;
	const int nx=2,ny=2,nz=2;
	int n[5]={ng,nv,nx,ny,nz};
	float *eta[ng*nv],*d_eta[nv];
	int gs[nx*ny*nz]={0};

	for(v1=0;v1<nv;v1++)
	{
		d_eta[v1]=new float[nx*ny*nz];
		for(g1=0;g1<ng;g1++)
		{		
			indx1=g1*nv+v1;
			eta[indx1]=new float[nx*ny*nz];
		}
	}

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=i*ny*nz+j*nz+k;
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=g1*nv+v1;
						(eta[indx1])[indx2]=0;
						(d_eta[v1])[indx2]=0;
						if(i<1 && j<1)
							(eta[gs[indx2]*nv+0])[indx2]=0;
					}
			}
	(eta[0])[0]=0.5;
	Ggradient aaa;
	aaa.Set(n,ka,ISO_P7);
	cout<<aaa.Energy(eta,gs)<<endl<<endl;
	aaa.Potential(eta,d_eta,gs);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=i*ny*nz+j*nz+k;
				for(v1=0;v1<nv;v1++)
				{
					for(g1=0;g1<ng;g1++)
					{
						indx1=g1*nv+v1;
						cout<<(eta[indx1])[indx2]<<'\t';
					}
				}
			}
	cout<<endl<<endl;
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=i*ny*nz+j*nz+k;
				for(v1=0;v1<nv;v1++)
					cout<<(d_eta[v1])[indx2]<<'\t';
			}
					
	for(v1=0;v1<nv;v1++)
	{
		delete[] d_eta[v1];
		for(g1=0;g1<ng;g1++)
		{		
			indx1=g1*nv+v1;
			delete[] eta[indx1];
		}
	}

	return 0;
}

int TryGInitConfig()
{
	int g1,v1;
	long int indx1;
	
	float a[1]={1};
	const int ng=1,nv=1;
	const int nx=3,ny=3,nz=3;
	int n[5]={ng,nv,nx,ny,nz};
	float *eta[ng*nv];
	int gs[nx*ny*nz];

	for(v1=0;v1<nv;v1++)
	{
		for(g1=0;g1<ng;g1++)
		{		
			indx1=g1*nv+v1;
			eta[indx1]=new float[nx*ny*nz];
		}
	}
	GInitConfig aaa;
	aaa.Set(n,SINGLE,eta,a,NULL,gs,NULL);

	return 0;
}

int TryGelastic()
{
	int n[5]={1,1,1,16,16};
	float c[3]={3.5,1.5,1};
	float a[3]={(float)1.0192,(float)0.8929,(float)1.0789};
	Gelastic aaa;
	aaa.Set(n,c,a,NULL,0,1,FIX,PHI1);
	aaa.Calc_Bpq();
	return 0;
}
*/
int TryGmt()
{
	int i,j,k;
	int nx,ny,nz;
/*
	ofstream fout("gs0.dat",ios::out);

	nx=128,ny=96,nz=72;
	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
			{
				if(i<nx/2)
					fout<<'0'<<endl;
				else
					fout<<'1'<<endl;

			}
	fout.flush();
	fout.close();
*/


	cout<<"Initializing..."<<endl;
	Gmt mt1("input.mt");
	cout<<"Evolution"<<endl;
         //mt1.Evolution();
//	mt1.Evolution_RK2();
        mt1.Evolution_Spectrum();
//	mt1.Output_eta();
//	mt1.Output_gs();
	return 0;
}

int gs_rand()
{
	int i,j,k,l;
	int nx,ny,nz;
	const int ng=6;
	nx=ny=nz=32;

	int x[ng],y[ng],z[ng];
	float d[ng]={0},tp1,tp2,dm=1000,gm=0;

	srand(time(0));
	for(l=0;l<ng;l++)
	{
		x[l]=rand()%nx;
		y[l]=rand()%ny;
		z[l]=rand()%nz;
	}
	
	ofstream fout("gs0.dat",ios::out);
	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
			{
				dm=1000;
				for(l=0;l<ng;l++)
				{
					d[l]=0;
					tp1=abs(i-x[l]);
					tp2=(tp1<nx/2)?tp1:nx-tp1;
					d[l]+=tp2*tp2;
					tp1=abs(j-y[l]);
					tp2=(tp1<ny/2)?tp1:ny-tp1;
					d[l]+=tp2*tp2;
					tp1=abs(k-z[l]);
					tp2=(tp1<nz/2)?tp1:nz-tp1;
					d[l]+=tp2*tp2;

					if(dm>d[l])
					{
						dm=d[l];
						gm=l;
					}
				}
				fout<<gm<<endl;
			}
	fout.flush();
	fout.close();


	return 0;
}

int main()
{
//	TryGfftw3();
//	TryGchemical();
//	TryGgradient();
//	TryGInitConfig();
//	TryGelastic();
	time_t time_begin,time_end;
	
	time(&time_begin);

	TryGmt();

	time(&time_end);
	cout<<ctime(&time_begin);
	cout<<ctime(&time_end)<<endl;

	ofstream tt("tt.txt",ios::out);
	tt<<ctime(&time_begin);
	tt<<ctime(&time_end)<<endl;
	tt.flush();
	tt.close();

	return 0;
}

