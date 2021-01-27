////////////////////////
// Gelastic.h	////////
// Elastic Energy	////
// Y.Gao 05-22-2009	////
////////////////////////

#include "Gelastic.h"
#include <iomanip>

Gelastic::Gelastic()
{
	int mm,nn;
	ng=nv=nx=ny=nz=0;
	c[0]=3.5;
	c[1]=1.5;
	c[2]=1;
	a[0]=a[1]=a[2]=1;
	g[0]=g[1]=g[2]=g_sqr=NULL;
	u=e0=s0=NULL;
	rot=NULL;
	Bpq=NULL;
	phi_k=NULL;
	s_app_miu=NULL;
	s_app_mag=0;
	e_app_miu=NULL;
	e_app_mag=0;
	scale=1;
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			Bpd[mm][nn]=NULL;
			s_app[mm][nn]=0;
			e_app[mm][nn]=0;
		}
	initflag=false;
	bc=FIX;
	ec=PHI1;
	cd=0;
	output_stress=0;
	output_interenergy=0;
}

Gelastic::~Gelastic()
{
	if(initflag==true)
		Destroy_Space();
}

int Gelastic::Init_Space()
{
	int g1,v1,g2,v2;
	int mm,nn;
	long int indx11,indx1,nxyz=nx*ny*nz;
//	Allocate e0,u,s0	
	e0=new float*[ng*nv];
	assert(e0!=NULL);
	u=new float*[ng*nv];
	assert(u!=NULL);
	s0=new float*[ng*nv];
	assert(s0!=NULL);
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			e0[indx1]=new float[9];
			assert(e0[indx1]!=NULL);
			u[indx1]=new float[9];
			assert(u[indx1]!=NULL);
			s0[indx1]=new float[9];
			assert(s0[indx1]!=NULL);
		}
//	Allocate Bpq,Bpd
//	if(output_stress!=1)
//	{
		Bpq=new float*[ng*nv*ng*nv];
		assert(Bpq!=NULL);
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
				for(g2=0;g2<ng;g2++)
					for(v2=0;v2<nv;v2++)
					{
						indx11=index11(g1,v1,g2,v2);
						Bpq[indx11]=new float[nxyz];
						assert(Bpq[indx11]!=NULL);
					}
/*	}
	else
	{
	*/
		for(mm=0;mm<3;mm++)
			for(nn=0;nn<3;nn++)
			{
				Bpd[mm][nn]=new float*[ng*nv];
				assert(Bpd[mm][nn]!=NULL);
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);	
						(Bpd[mm][nn])[indx1]=new float[nx*ny*nz];
						assert((Bpd[mm][nn])[indx1]!=NULL);
					}
			}
//	}
//	Allocate rot
	rot=new float*[ng];
	assert(rot!=NULL);
	for(g1=0;g1<ng;g1++)
	{
		rot[g1]=new float[9];
		assert(rot[g1]!=NULL);
	}
//	Allocate g,g_sqr
	g_sqr=new float[nxyz];
	assert(g_sqr!=NULL);
	for(mm=0;mm<3;mm++)
	{
		g[mm]=new float[nxyz];
		assert(g[mm]!=NULL);
	}
//	Allocate s_app_miu,e_app_miu
	s_app_miu=new float[ng*nv];
	assert(s_app_miu!=NULL);
	e_app_miu=new float[ng*nv];
	assert(e_app_miu!=NULL);
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			s_app_miu[indx1]=0;
			e_app_miu[indx1]=0;
		}
//	Allocate stress field
	s_field=new float*[9];
	assert(s_field!=NULL);
	for(int i=0;i<9;i++)
	{
		s_field[i]=new float[nxyz];
		assert(s_field[i]!=NULL);
	}
//	Allocate phi_k,miu
	phi_k=new fftw_complex*[ng*nv];
	assert(phi_k!=NULL);
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			
			phi_k[indx1]=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nxyz);
			assert(phi_k[indx1]!=NULL);
		}

	return 0;
}

int Gelastic::Destroy_Space()
{
	int g1,v1;
	int mm;
	long int indx1;
//	Deallocate e0,u,s0
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			delete[] e0[indx1];
			delete[] u[indx1];
			delete[] s0[indx1];
		}
	delete[] e0;
	delete[] u;
	delete[] s0;
	e0=NULL;
	u=NULL;
	s0=NULL;
//	Deallocate Bpq
	Destroy_Bpq();
//  Deallocate rot
	for(g1=0;g1<ng;g1++)
		delete[] rot[g1];
	delete[] rot;
	rot=NULL;
//	Deallocate g,g_sqr
	delete[] g_sqr;
	g_sqr=NULL;
	for(mm=0;mm<3;mm++)
	{
		delete[] g[mm];
		g[mm]=NULL;
	}
//	Deallocate s_app_miu
	delete[] s_app_miu;
	s_app_miu=NULL;
//	Deallocate e_app_miu
	delete[] e_app_miu;
	e_app_miu=NULL;
//	Deallocate phi_k
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			fftw_free(phi_k[indx1]);
		}
	delete[] phi_k;
	phi_k=NULL;

	return 0;
}

int Gelastic::Destroy_Bpq()
{
	int g1,v1,g2,v2;
	long int indx11,indx1;
	int mm,nn;

	if(output_stress!=1)
	{
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
				for(g2=0;g2<ng;g2++)
					for(v2=0;v2<nv;v2++)
					{
						indx11=index11(g1,v1,g2,v2);
						delete[] Bpq[indx11];
					}
		delete[] Bpq;
		Bpq=NULL;
	}
	else
	{
		for(mm=0;mm<3;mm++)
			for(nn=0;nn<3;nn++)
			{
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);	
						delete[] (Bpd[mm][nn])[indx1];
					}
				delete[] Bpd[mm][nn];
				Bpd[mm][nn]=NULL;
			}
	}

	return 0;
}

int Gelastic::Set(int n[],float cc[],float aa[],
				  float **rott,float s_app_magg,float s_app_incc, float e_app_magg,float scalee,
				  boundary_condition bcc,elastic_couple ecc,int cdd,fftw_option op,
				  int output_s,int output_ie)
{
	s_app_mag=s_app_magg;
	e_app_mag=e_app_magg;
	s_app_inc=s_app_incc;
	scale=scalee;
	bc=bcc;
	ec=ecc;
	cd=cdd;
	output_stress=output_s;
	output_interenergy=output_ie;
	initflag=true;
	
	Set_n(n);
	Init_Space();

	gfft.Set(n,op);

	Set_g();
	Set_Cijkl(cc);
	Set_a(aa);
	Set_SFTS(); // w/o rotation
	Set_rot(rott);
	Set_s0();
        //Calc_Bpd();

	if(output_interenergy==1)
		return 1;
	if(output_stress!=1)
{		Calc_Bpd();
		Calc_Bpq();
	}else
		Calc_Bpd();

	return 0;
}

int Gelastic::Set_n(int n[])
{
	if(n[0]<1 || n[1]<1 || n[2]<1 || n[3]<1 || n[4]<1)
	{
		cout<<"Invalid dimensions in Gelastic!"<<endl;
		exit(1);
	}
	ng=n[0];
	nv=n[1];
	nx=n[2];
	ny=n[3];
	nz=n[4];

	return 0;
}

int Gelastic::Set_a(float aa[])
{
	int i;
	for(i=0;i<4;i++)
		a[i]=aa[i];
	return 0;
}

int Gelastic::Set_Cijkl(float cc[])
{
	int ii,jj,kk,ll;
	
	c[0]=cc[0];//c11
	c[1]=cc[1];//c12
	c[2]=cc[2];//c44
	
	for(ii=0;ii<3;ii++)
		for(jj=0;jj<3;jj++)
			for(kk=0;kk<3;kk++)
				for(ll=0;ll<3;ll++)
				{
					Cijkl[ii][jj][kk][ll]=0;
				}
	Cijkl[0][0][0][0]=Cijkl[1][1][1][1]=Cijkl[2][2][2][2]
		=c[0];
	Cijkl[0][0][1][1]=Cijkl[1][1][0][0]
		=Cijkl[1][1][2][2]=Cijkl[2][2][1][1]
		=Cijkl[2][2][0][0]=Cijkl[0][0][2][2]
		=c[1];
	Cijkl[0][1][0][1]=Cijkl[0][1][1][0]=Cijkl[1][0][0][1]=Cijkl[1][0][1][0]
		=Cijkl[1][2][1][2]=Cijkl[1][2][2][1]=Cijkl[2][1][1][2]=Cijkl[2][1][2][1]
		=Cijkl[2][0][2][0]=Cijkl[2][0][0][2]=Cijkl[0][2][2][0]=Cijkl[0][2][0][2]
		=c[2];

//	Calculate Sijkl

	return 0;
}

int Gelastic::Set_SFTS()
{
	int g1,v1;
	long int indx1;
	int ii,jj,kk,iindx2;

	float eye[3][3]={0};
	eye[0][0]=eye[1][1]=eye[2][2]=1;
	
	Set_u();

//	SFTS w/o rotation	
	/*
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				{
					iindx2=ii*3+jj;
					(e0[indx1])[iindx2]=0;
					for(kk=0;kk<3;kk++)
						(e0[indx1])[iindx2]+=(u[indx1])[kk*3+ii]*(u[indx1])[kk*3+jj];

					(e0[indx1])[iindx2]-=eye[ii][jj];
					(e0[indx1])[iindx2]/=2.0;
				}
		}
*/
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
//               if(v1=12)
  //               continue;
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				{
					iindx2=ii*3+jj;
					(e0[indx1])[iindx2]=0;
              // if(v1>3)
               //  continue;
					for(kk=0;kk<3;kk++)
						(e0[indx1])[iindx2]+=(u[indx1])[kk*3+ii]*(u[indx1])[kk*3+jj];

					(e0[indx1])[iindx2]-=eye[ii][jj];
					(e0[indx1])[iindx2]/=2.0;
                        }
		}
//	SFTS_habit for sivom
//	Set_SFTS_habit();

//  SFTS for NiPtTi prep
//	Set_SFTS_NiPtTi();

//	SFTS for Laminate structure of P phase
//	Set_SFTS_PLaminate();

//	SFTS for TiAlV
//	Set_SFTS_TiAlV();

//	SFTS for H-phase
//	Set_SFTS_NiTiHf();
//	set SFTS for ppt
// use B2 to tetragonal
/*
        e0[0][0]=-0.0906;
        e0[0][4]=-0.0906;
        e0[0][5]=0.2895;
        e0[0][7]=0.2895;
        e0[0][8]=-0.0906;
  */     
 
//        e0[0][0]=0.1;
//        e0[0][4]=-0.05;
//        e0[0][8]=-0.05;
//        //e0[0][0]=2;
//        e0[1][4]=0.1;
//        e0[1][0]=-0.05;
//        e0[1][8]=-0.05;
//        
//        e0[2][8]=0.1;
//        e0[2][0]=-0.05;
//        e0[2][4]=-0.05;

/*
        e0[12][0]=-0.1;
        e0[12][4]=-0.1;
        e0[12][5]=0.0;
        e0[12][7]=0.0;
        e0[12][8]=-0.1;
*/
//        e0[0][0]=0.1;
//        e0[0][4]=-0.1;
//
//        e0[1][0]=-0.1;
//        e0[1][4]=0.1;

/*
				e0[12][0]=0.00195;
				e0[12][1]=0.00726;
				e0[12][2]=0;
				e0[12][3]=0.00726;
				e0[12][4]=0.00195;
				e0[12][5]=0.0;
				e0[12][6]=0;
				e0[12][7]=0.0;
				e0[12][8]=0.02291;
 */
//SFTS of ppt
//
//
/*
                                e0[4][0]=0.0224;
				e0[4][1]=0;
				e0[4][2]=0;
				e0[4][3]=0;
				e0[4][4]=0.0024;
				e0[4][5]=0.0072;
				e0[4][6]=0;
				e0[4][7]=0.0072;
				e0[4][8]=0.0024;


				e0[5][0]=0.0224;
				e0[5][1]=0;
				e0[5][2]=0;
				e0[5][3]=0;
				e0[5][4]=0.0024;
				e0[5][5]=-0.0072;
				e0[5][6]=0;
				e0[5][7]=-0.0072;
				e0[5][8]=0.0024;
				


				e0[6][0]=0.0024;
				e0[6][1]=0;
				e0[6][2]=0.0072;
				e0[6][3]=0;
				e0[6][4]=0.0224;
				e0[6][5]=0;
				e0[6][6]=0.0072;
				e0[6][7]=0;
				e0[6][8]=0.0024;

				e0[7][0]=0.0024;
				e0[7][1]=0;
				e0[7][2]=-0.0072;
				e0[7][3]=0;
				e0[7][4]=0.0224;
				e0[7][5]=0;
				e0[7][6]=-0.0072;
				e0[7][7]=0;
				e0[7][8]=0.0024;
				
				e0[8][0]=0.0024;
				e0[8][1]=0.0072;
				e0[8][2]=0;
				e0[8][3]=0.0072;
				e0[8][4]=0.0024;
				e0[8][5]=0;
				e0[8][6]=0;
				e0[8][7]=0;
				e0[8][8]=0.0224;
				
				e0[9][0]=0.0024;
				e0[9][1]=-0.0072;
				e0[9][2]=0;
				e0[9][3]=-0.0072;
				e0[9][4]=0.0024;
				e0[9][5]=0;
				e0[9][6]=0;
				e0[9][7]=0;
				e0[9][8]=0.0224;
*/

/*
				e0[4][0]=0.0;
				e0[4][1]=0;
				e0[4][2]=0;
				e0[4][3]=0;
				e0[4][4]=0.0;
				e0[4][5]=0.0;
				e0[4][6]=0;
				e0[4][7]=0.0;
				e0[4][8]=0.0;


				e0[5][0]=0.0;
				e0[5][1]=0;
				e0[5][2]=0;
				e0[5][3]=0;
				e0[5][4]=0.0;
				e0[5][5]=-0.0;
				e0[5][6]=0;
				e0[5][7]=-0.0;
				e0[5][8]=0.0;
				


				e0[6][0]=0.0;
				e0[6][1]=0;
				e0[6][2]=0.0;
				e0[6][3]=0;
				e0[6][4]=0.0;
				e0[6][5]=0;
				e0[6][6]=0.0;
				e0[6][7]=0;
				e0[6][8]=0.0;

				e0[7][0]=0.0;
				e0[7][1]=0;
				e0[7][2]=-0.0;
				e0[7][3]=0;
				e0[7][4]=0.0;
				e0[7][5]=0;
				e0[7][6]=-0.0;
				e0[7][7]=0;
				e0[7][8]=0.0;
				
				e0[8][0]=0.0;
				e0[8][1]=0.0;
				e0[8][2]=0;
				e0[8][3]=0.0;
				e0[8][4]=0.0;
				e0[8][5]=0;
				e0[8][6]=0;
				e0[8][7]=0;
				e0[8][8]=0.0;
				
				e0[9][0]=0.0;
				e0[9][1]=-0.0;
				e0[9][2]=0;
				e0[9][3]=-0.0;
				e0[9][4]=0.0;
				e0[9][5]=0;
				e0[9][6]=0;
				e0[9][7]=0;
				e0[9][8]=0.0;

*/
/* 
				e0[10][0]=0.0;
				e0[10][1]=-0.0;
				e0[10][2]=0;
				e0[10][3]=-0.00;
				e0[10][4]=0.00;
				e0[10][5]=0;
				e0[10][6]=0;
				e0[10][7]=0;
				e0[10][8]=0.0;

				e0[11][0]=0.0;
				e0[11][1]=-0.0;
				e0[11][2]=0;
				e0[11][3]=-0.00;
				e0[11][4]=0.00;
				e0[11][5]=0;
				e0[11][6]=0;
				e0[11][7]=0;
				e0[11][8]=0.0;
*/
return 0;
}

int Gelastic::Set_u()
{
	int g1,v1;
	long int indx1;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			switch(v1)
			{
			case 0:
				// Orthorhombic
				/*
				(u[indx1])[0]=(a[0]+a[2])/2.0;
				(u[indx1])[1]=0;
				(u[indx1])[2]=(a[0]-a[2])/2.0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=0;
				(u[indx1])[6]=(a[0]-a[2])/2.0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=(a[0]+a[2])/2.0;
				*/
/*
				// Tetragonal
				(u[indx1])[0]=sqrt(1.4);
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=sqrt(0.8);
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=sqrt(0.8);

				(u[indx1])[0]=1.001;
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=0.999;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=1;
*/
				// Monoclinic

				(u[indx1])[0]=a[0];
				(u[indx1])[1]=a[2];
				(u[indx1])[2]=a[2];
				(u[indx1])[3]=a[2];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=a[3];
				(u[indx1])[6]=a[2];
				(u[indx1])[7]=a[3];
				(u[indx1])[8]=a[1];
				
				// Averaged Trigonal
/*				(u[indx1])[0]=0.9842;
				(u[indx1])[1]=-0.0017;
				(u[indx1])[2]=-0.0017;
				(u[indx1])[3]=-0.0017;
				(u[indx1])[4]=0.9842;
				(u[indx1])[5]=-0.0017;
				(u[indx1])[6]=-0.0017;
				(u[indx1])[7]=-0.0017;
				(u[indx1])[8]=0.9842;
*/

				break;
			case 1:
				// Orthorhombic
/*				(u[indx1])[0]=(a[0]+a[2])/2.0;
				(u[indx1])[1]=0;
				(u[indx1])[2]=-(a[0]-a[2])/2.0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=0;
				(u[indx1])[6]=-(a[0]-a[2])/2.0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=(a[0]+a[2])/2.0;
*/
				// Tetragonal
/*				(u[indx1])[0]=sqrt(0.8);
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=sqrt(1.4);
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=sqrt(0.8);

				(u[indx1])[0]=0.8;
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=1.5;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=1;				
*/
				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=a[2];
				(u[indx1])[2]=a[3];
				(u[indx1])[3]=a[2];
				(u[indx1])[4]=a[0];
				(u[indx1])[5]=a[2];
				(u[indx1])[6]=a[3];
				(u[indx1])[7]=a[2];
				(u[indx1])[8]=a[1];				

				// Averaged Trigonal
/*				(u[indx1])[0]=0.9842;
				(u[indx1])[1]=0.0017;
				(u[indx1])[2]=0.0017;
				(u[indx1])[3]=0.0017;
				(u[indx1])[4]=0.9842;
				(u[indx1])[5]=-0.0017;
				(u[indx1])[6]=0.0017;
				(u[indx1])[7]=-0.0017;
				(u[indx1])[8]=0.9842;
*/
				break;
			case 2:
				// Orthorhombic
/*				(u[indx1])[0]=(a[0]+a[2])/2.0;
				(u[indx1])[1]=(a[0]-a[2])/2.0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=(a[0]-a[2])/2.0;
				(u[indx1])[4]=(a[0]+a[2])/2.0;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=a[1];
*/
/*
				// Tetragonal
				(u[indx1])[0]=0.913;
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=0.913;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=1.2;
*/
				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=a[3];
				(u[indx1])[2]=a[2];
				(u[indx1])[3]=a[3];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=a[2];
				(u[indx1])[6]=a[2];
				(u[indx1])[7]=a[2];
				(u[indx1])[8]=a[0];				

				// Averaged Trigonal
/*				(u[indx1])[0]=0.9842;
				(u[indx1])[1]=0.0017;
				(u[indx1])[2]=-0.0017;
				(u[indx1])[3]=0.0017;
				(u[indx1])[4]=0.9842;
				(u[indx1])[5]=0.0017;
				(u[indx1])[6]=-0.0017;
				(u[indx1])[7]=0.0017;
				(u[indx1])[8]=0.9842;
*/
				break;

			case 3:
				// Orthorhombic
/*				(u[indx1])[0]=(a[0]+a[2])/2.0;
				(u[indx1])[1]=-(a[0]-a[2])/2.0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=-(a[0]-a[2])/2.0;
				(u[indx1])[4]=(a[0]+a[2])/2.0;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=a[1];
*/

				// Monoclinic
				(u[indx1])[0]=a[0];
				(u[indx1])[1]=-a[2];
				(u[indx1])[2]=-a[2];
				(u[indx1])[3]=-a[2];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=a[3];
				(u[indx1])[6]=-a[2];
				(u[indx1])[7]=a[3];
				(u[indx1])[8]=a[1];

				// Averaged Trigonal
			/*	(u[indx1])[0]=0.9842;
				(u[indx1])[1]=-0.0017;
				(u[indx1])[2]=0.0017;
				(u[indx1])[3]=-0.0017;
				(u[indx1])[4]=0.9842;
				(u[indx1])[5]=0.0017;
				(u[indx1])[6]=0.0017;
				(u[indx1])[7]=0.0017;
				(u[indx1])[8]=0.9842;
*/
				break;
			case 4:
				// Orthorhombic
				/*
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=(a[0]+a[2])/2.0;
				(u[indx1])[5]=(a[0]-a[2])/2.0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=(a[0]-a[2])/2.0;
				(u[indx1])[8]=(a[0]+a[2])/2.0;
*/

				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=-a[2];
				(u[indx1])[2]=-a[3];
				(u[indx1])[3]=-a[2];
				(u[indx1])[4]=a[0];
				(u[indx1])[5]=a[2];
				(u[indx1])[6]=-a[3];
				(u[indx1])[7]=a[2];
				(u[indx1])[8]=a[1];				
		
				break;
			case 5:
				// Orthorhombic
				/*
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=(a[0]+a[2])/2.0;
				(u[indx1])[5]=-(a[0]-a[2])/2.0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=-(a[0]-a[2])/2.0;
				(u[indx1])[8]=(a[0]+a[2])/2.0;
*/

				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=-a[3];
				(u[indx1])[2]=-a[2];
				(u[indx1])[3]=-a[3];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=a[2];
				(u[indx1])[6]=-a[2];
				(u[indx1])[7]=a[2];
				(u[indx1])[8]=a[0];				
			
				break;
			case 6:
				// Monoclinic
				(u[indx1])[0]=a[0];
				(u[indx1])[1]=-a[2];
				(u[indx1])[2]=a[2];
				(u[indx1])[3]=-a[2];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=-a[3];
				(u[indx1])[6]=a[2];
				(u[indx1])[7]=-a[3];
				(u[indx1])[8]=a[1];				
				break;
			case 7:
				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=-a[2];
				(u[indx1])[2]=a[3];
				(u[indx1])[3]=-a[2];
				(u[indx1])[4]=a[0];
				(u[indx1])[5]=-a[2];
				(u[indx1])[6]=a[3];
				(u[indx1])[7]=-a[2];
				(u[indx1])[8]=a[1];	
				break;
			case 8:
				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=-a[3];
				(u[indx1])[2]=a[2];
				(u[indx1])[3]=-a[3];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=-a[2];
				(u[indx1])[6]=a[2];
				(u[indx1])[7]=-a[2];
				(u[indx1])[8]=a[0];				
				break;
			case 9:
				// Monoclinic
				(u[indx1])[0]=a[0];
				(u[indx1])[1]=a[2];
				(u[indx1])[2]=-a[2];
				(u[indx1])[3]=a[2];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=-a[3];
				(u[indx1])[6]=-a[2];
				(u[indx1])[7]=-a[3];
				(u[indx1])[8]=a[1];
				break;
			case 10:
				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=a[2];
				(u[indx1])[2]=-a[3];
				(u[indx1])[3]=a[2];
				(u[indx1])[4]=a[0];
				(u[indx1])[5]=-a[2];
				(u[indx1])[6]=-a[3];
				(u[indx1])[7]=-a[2];
				(u[indx1])[8]=a[1];				
				break;
			case 11:
				// Monoclinic
				(u[indx1])[0]=a[1];
				(u[indx1])[1]=a[3];
				(u[indx1])[2]=-a[2];
				(u[indx1])[3]=a[3];
				(u[indx1])[4]=a[1];
				(u[indx1])[5]=-a[2];
				(u[indx1])[6]=-a[2];
				(u[indx1])[7]=-a[2];
				(u[indx1])[8]=a[0];				
				break;
/*			case 12:
				// Concentration
				(u[indx1])[0]=1;
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=1;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=1;				
				break;
*/
			default:
				(u[indx1])[0]=1;
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=1;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=1;				
				break;
			}
		}
	return 0;
}

int Gelastic::Set_s0()
{
	int g1,v1;
	long int indx1;
	int ii,jj,kk,ll,iindx2;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				{
					iindx2=ii*3+jj;
					(s0[indx1])[iindx2]=0;
					for(kk=0;kk<3;kk++)
						for(ll=0;ll<3;ll++)
							(s0[indx1])[iindx2]+=
							Cijkl[ii][jj][kk][ll]*(e0[indx1])[kk*3+ll];
				}
			
		}
	return 0;
}

int Gelastic::Set_rot(float **rott)
{
	int g1,v1,ii,jj,mm,nn,indx2;
	long int indx1;

	// No rotation?
	g1=0;
	if(g1<ng)
	{

//		(111)->(100)
/*
		(rot[g1])[0]=sqrt(3.0)/3.0;
		(rot[g1])[1]=sqrt(3.0)/3.0;
		(rot[g1])[2]=sqrt(3.0)/3.0;
		(rot[g1])[3]=sqrt(2.0)/2.0;
		(rot[g1])[4]=-sqrt(2.0)/2.0;
		(rot[g1])[5]=0;
		(rot[g1])[6]=sqrt(6.0)/6.0;
		(rot[g1])[7]=sqrt(6.0)/6.0;
		(rot[g1])[8]=-sqrt(6.0)/3.0;
*/
/*
		// z-45
		(rot[g1])[0]=sqrt(2)/2.0;
		(rot[g1])[1]=sqrt(2)/2.0;
		(rot[g1])[2]=0;
		(rot[g1])[3]=-sqrt(2)/2.0;
		(rot[g1])[4]=sqrt(2)/2.0;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;
*/
		// no rotation
	
		(rot[g1])[0]=1;
		(rot[g1])[1]=0;
		(rot[g1])[2]=0;
		(rot[g1])[3]=0;
		(rot[g1])[4]=1;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;


		// rotaiton for H-phase {113} facet
/*
		(rot[g1])[0]=0.7071;
		(rot[g1])[1]=0;
		(rot[g1])[2]=0.7071;
		(rot[g1])[3]=0;
		(rot[g1])[4]=1;
		(rot[g1])[5]=0;
		(rot[g1])[6]=-0.7071;
		(rot[g1])[7]=0;
		(rot[g1])[8]=0.7071;
*/
/*
		(rot[g1])[0]=0.0;
		(rot[g1])[1]=-0.7071;
		(rot[g1])[2]=0.7071;
		(rot[g1])[3]=1;
		(rot[g1])[4]=0;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0.0;
		(rot[g1])[7]=0.7071;
		(rot[g1])[8]=0.7071;
*/
		cout<<"grain 0 is rotated!"<<endl;
}

	// rotate 10 degree along z axis
	g1=1;
	if(g1<ng)
	{
		(rot[g1])[0]=(float)0.9848; // z
		(rot[g1])[1]=(float)0.1736;
		(rot[g1])[2]=0;
		(rot[g1])[3]=(float)-0.1736;
		(rot[g1])[4]=(float)0.9848;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;

		(rot[g1])[0]=(float)0.9848; // y
		(rot[g1])[1]=0;
		(rot[g1])[2]=(float)0.1736;
		(rot[g1])[3]=0;
		(rot[g1])[4]=1;
		(rot[g1])[5]=0;
		(rot[g1])[6]=(float)-0.1736;
		(rot[g1])[7]=0;
		(rot[g1])[8]=(float)0.9848;
/*
		(rot[g1])[0]=1; // x
		(rot[g1])[1]=0;
		(rot[g1])[2]=0;
		(rot[g1])[3]=0;
		(rot[g1])[4]=(float)0.9848;
		(rot[g1])[5]=(float)0.1736;
		(rot[g1])[6]=0;
		(rot[g1])[7]=(float)-0.1736;
		(rot[g1])[8]=(float)0.9848;
*/	}

	// rotate 20 degree along z axis
	g1=2;
	if(g1<ng)
	{
		(rot[g1])[0]=(float)0.9397;
		(rot[g1])[1]=(float)0.3420;
		(rot[g1])[2]=0;
		(rot[g1])[3]=(float)-0.3420;
		(rot[g1])[4]=(float)0.9397;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;
	}

	// rotate 30 degree along z axis
	g1=3;
	if(g1<ng)
	{
		(rot[g1])[0]=(float)0.8660;
		(rot[g1])[1]=(float)0.5;
		(rot[g1])[2]=0;
		(rot[g1])[3]=(float)-0.5;
		(rot[g1])[4]=(float)0.866;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;
	}

	// rotate 40 degree along z axis
	g1=4;
	if(g1<ng)
	{
		(rot[g1])[0]=(float)0.7660;
		(rot[g1])[1]=(float)0.6428;
		(rot[g1])[2]=0;
		(rot[g1])[3]=(float)-0.6428;
		(rot[g1])[4]=(float)0.7660;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;

		(rot[g1])[0]=(float)0.7660; // y
		(rot[g1])[1]=0;
		(rot[g1])[2]=(float)0.6428;
		(rot[g1])[3]=0;
		(rot[g1])[4]=1;
		(rot[g1])[5]=0;
		(rot[g1])[6]=(float)-0.6428;
		(rot[g1])[7]=0;
		(rot[g1])[8]=(float)0.7660;

		(rot[g1])[0]=1; // x
		(rot[g1])[1]=0;
		(rot[g1])[2]=0;
		(rot[g1])[3]=0;
		(rot[g1])[4]=(float)0.7660;
		(rot[g1])[5]=(float)0.6428;
		(rot[g1])[6]=0;
		(rot[g1])[7]=(float)-0.6428;
		(rot[g1])[8]=(float)0.7660;
}

	// rotate 45 degree along z axis
	g1=5;
	if(g1<ng)
	{
		(rot[g1])[0]=(float)0.707;
		(rot[g1])[1]=(float)0.707;
		(rot[g1])[2]=0;
		(rot[g1])[3]=(float)-0.707;
		(rot[g1])[4]=(float)0.707;
		(rot[g1])[5]=0;
		(rot[g1])[6]=0;
		(rot[g1])[7]=0;
		(rot[g1])[8]=1;
	}
	// Keep grain 0, rotate other grain
	for(g1=1;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				{
					indx2=ii*3+jj;
					(e0[indx1])[indx2]=0;
					for(mm=0;mm<3;mm++)
						for(nn=0;nn<3;nn++)
							(e0[indx1])[indx2]+=(rot[g1])[ii*3+mm]*(rot[g1])[jj*3+nn]*(e0[v1])[mm*3+nn];
				}
		}

// Rotate stain and Cijkl grain 0
	g1=0;
	float e_tp[3][3]={0},Cijkl_tp[3][3][3][3]={0};
	//strain
	for(v1=0;v1<nv;v1++)
	{
		indx1=index1(g1,v1);
		for(ii=0;ii<3;ii++)
			for(jj=0;jj<3;jj++)
			{
				e_tp[ii][jj]=0;
				for(mm=0;mm<3;mm++)
					for(nn=0;nn<3;nn++)
						e_tp[ii][jj]+=rot[g1][ii*3+mm]*rot[g1][jj*3+nn]*e0[indx1][mm*3+nn];
			}
		for(ii=0;ii<3;ii++)
			for(jj=0;jj<3;jj++)
			{
				e0[indx1][ii*3+jj]=e_tp[ii][jj];
			}
	}
	int ip,jp,kp,lp,kk,ll;
	for(ip=0;ip<3;ip++)
	for(jp=0;jp<3;jp++)
	for(kp=0;kp<3;kp++)
	for(lp=0;lp<3;lp++)
	{
		Cijkl_tp[ip][jp][kp][lp]=0;
		for(ii=0;ii<3;ii++)
		for(jj=0;jj<3;jj++)
		for(kk=0;kk<3;kk++)
		for(ll=0;ll<3;ll++)
		{
			Cijkl_tp[ip][jp][kp][lp]+=rot[g1][ip*3+ii]*rot[g1][jp*3+jj]*rot[g1][kp*3+kk]*rot[g1][lp*3+ll]*Cijkl[ii][jj][kk][ll];
		}
	}
	for(ip=0;ip<3;ip++)
	for(jp=0;jp<3;jp++)
	for(kp=0;kp<3;kp++)
	for(lp=0;lp<3;lp++)
	{
		Cijkl[ip][jp][kp][lp]=Cijkl_tp[ip][jp][kp][lp];
		cout<<ip+1<<' '<<jp+1<<' '<<kp+1<<' '<<lp+1<<'\t'<<Cijkl[ip][jp][kp][lp]<<endl;
	}

//	For H-phase test
/*
	e0[0][0]=0.00195;
	e0[0][1]=0.00726;
	e0[0][2]=0;
	e0[0][3]=0.00726;
	e0[0][4]=0.00195;
	e0[0][5]=0;
	e0[0][6]=0;
	e0[0][7]=0;
	e0[0][8]=0.02291;
	e0[1][0]=0.00195;
	e0[1][1]=-0.00726;
	e0[1][2]=0;
	e0[1][3]=-0.00726;
	e0[1][4]=0.00195;
	e0[1][5]=0;
	e0[1][6]=0;
	e0[1][7]=0;
	e0[1][8]=0.02291;
	e0[2][0]=0.00195;
	e0[2][1]=0.0;
	e0[2][2]=0.00726;
	e0[2][3]=0.0;
	e0[2][4]=0.02291;
	e0[2][5]=0;
	e0[2][6]=0.00726;
	e0[2][7]=0;
	e0[2][8]=0.00195;
	e0[3][0]=0.00195;
	e0[3][1]=0.0;
	e0[3][2]=-0.00726;
	e0[3][3]=0.0;
	e0[3][4]=0.02291;
	e0[3][5]=0;
	e0[3][6]=-0.00726;
	e0[3][7]=0;
	e0[3][8]=0.00195;
	e0[4][0]=0.02291;
	e0[4][1]=0.0;
	e0[4][2]=0.0;
	e0[4][3]=0.0;
	e0[4][4]=0.00195;
	e0[4][5]=0.00726;
	e0[4][6]=0.0;
	e0[4][7]=0.00726;
	e0[4][8]=0.00195;
	e0[5][0]=0.02291;
	e0[5][1]=0.0;
	e0[5][2]=0.0;
	e0[5][3]=0.0;
	e0[5][4]=0.00195;
	e0[5][5]=-0.00726;
	e0[5][6]=0.0;
	e0[5][7]=-0.00726;
	e0[5][8]=0.00195;
*/
// Output SFTS after rotation
	for(g1=0;g1<ng;g1++)
	{
		cout<<"In grain "<<g1<<endl;
		for(v1=0;v1<nv;v1++)
		{
			cout<<"variant "<<v1<<endl;
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
			{
				for(jj=0;jj<3;jj++)
					cout<<e0[indx1][ii*3+jj]<<'\t';
				cout<<endl;
			}
			cout<<endl;
		}
		cout<<endl;
	}

	return 0;
}

float Gelastic::Energy()
{
	float result=0;
	return result;
}

int Gelastic::Potential(float **eta,float **d_eta,int *gs)
{
	int g1,v1,g2,v2,i,j,k;
	long int indx1_1,indx1_2,indx2,indx11;
	float e1,e2,e3,e4;
	float tp;
	
	if(fabs(scale)<1E-5)
		return 1;
	
	Calc_phi_k(eta);

	fftw_complex *miu_k;
	miu_k=new fftw_complex[nx*ny*nz];
	assert(miu_k!=NULL);
	float *miu=new float[nx*ny*nz];
	assert(miu!=NULL);

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1_1=index1(g1,v1);
			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);

						miu_k[indx2][0]=0;
						miu_k[indx2][1]=0;

						for(g2=0;g2<ng;g2++)
							for(v2=0;v2<nv;v2++)
							{
								indx1_2=index1(g2,v2);
								indx11=index11(g1,v1,g2,v2);

								miu_k[indx2][0]+=(Bpq[indx11])[indx2]*(phi_k[indx1_2])[indx2][0];
								miu_k[indx2][1]+=(Bpq[indx11])[indx2]*(phi_k[indx1_2])[indx2][1];
							}
					}
			gfft.FFTW_3D(miu_k,miu,BACKWARD);
			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);
						if(g1!=gs[indx2])
							continue;
						e1=eta[indx1_1][indx2];
						e2=e1*e1;
						e3=e2*e1;
						e4=e3*e1;
						switch(ec)
						{
						case PHI1:
							(d_eta[v1])[indx2]+=miu[indx2]*scale;
							break;
						case PHI2:
							(d_eta[v1])[indx2]+=2.0*e1*miu[indx2]*scale;
							break;
						case PHI23:
							(d_eta[v1])[indx2]+=6*(-e2+e1)*miu[indx2]*scale;
							break;
						case PHI345:
//							if(v1<4 || v1>5)
								d_eta[v1][indx2]+=30*(e4-2*e3+e2)*miu[indx2]*scale;
// Used for multi-field
/*							else if(v1==4)
							{
								tp=cosh((e1-0.3)*26.667);
								tp=1.0/tp/tp;
								d_eta[v1][indx2]+=tp*miu[indx2]*scale;
							}
							else
							{
								tp=cosh((e1-0.2)*60.606);
								tp=1.0/tp/tp;
								d_eta[v1][indx2]+=tp*miu[indx2]*scale;
							}
*/
							break;
						default:
							cout<<"Invalid elastic couple!"<<endl;
							exit(1);
						}
						if(fabs(s_app_mag)>1E-5)
//							(d_eta[v1])[indx2]+=s_app_miu[indx1_1]*6*(-e2+e1);
							(d_eta[v1])[indx2]+=s_app_miu[indx1_1]*2*(e1);
						if(fabs(e_app_mag)>1E-5)
//							(d_eta[v1])[indx2]+=e_app_miu[indx1_1]*6*(-e2+e1);
							(d_eta[v1])[indx2]+=e_app_miu[indx1_1]*2*(e1);
					}			
		}

	delete[] miu_k;
	delete[] miu;

	return 0;
}

int Gelastic::Potential_k(fftw_complex **eta,fftw_complex **d_eta_k) // used only elastic couple= linear
{
	if(fabs(s_app_mag)>1E-5)
	{
		cout<<"Invalid applied stress option in Potential_k!"<<endl;
		exit(1);
	}
	if(fabs(e_app_mag)>1E-5)
	{
		cout<<"Invalid applied strain option in Potential_k!"<<endl;
		exit(1);
	}

	
	return 0;
}

int Gelastic::Calc_phi_k(float **eta)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	float *temp;
	float e1,e2,e3,e4,e5;

	temp=new float[nx*ny*nz];

	switch(ec)
	{
	case PHI1:
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);
				gfft.FFTW_3D(eta[indx1],phi_k[indx1],FORWARD);
			}
		break;
	case PHI2:
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);
				for(i=0;i<nx;i++)
					for(j=0;j<ny;j++)
						for(k=0;k<nz;k++)
						{
							indx2=index2(i,j,k);
							e1=(eta[indx1])[indx2];
							e2=e1*e1;
							temp[indx2]=e2;
						}
				gfft.FFTW_3D(temp,phi_k[indx1],FORWARD);
			}
		break;
	case PHI23:
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);
				for(i=0;i<nx;i++)
					for(j=0;j<ny;j++)
						for(k=0;k<nz;k++)
						{
							indx2=index2(i,j,k);
							e1=(eta[indx1])[indx2];
							e2=e1*e1;
							e3=e2*e1;
							temp[indx2]=-2*e3+3*e2;
						}
				gfft.FFTW_3D(temp,phi_k[indx1],FORWARD);
			}
		break;
	case PHI345:
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);
				for(i=0;i<nx;i++)
					for(j=0;j<ny;j++)
						for(k=0;k<nz;k++)
						{
							indx2=index2(i,j,k);
							e1=eta[indx1][indx2];
							e2=e1*e1;
							e3=e2*e1;
							e4=e3*e1;
							e5=e4*e1;
//							if(v1<4 || v1>5)
								temp[indx2]=6*e5-15*e4+10*e3;
// Used for multi-field
/*							else if(v1==4)
								temp[indx2]=tanh((e1-0.3)*26.667)/26.667;
							else if(v1==5)
								temp[indx2]=tanh((e1-0.2)*60.606)/60.606;
*/						}
				gfft.FFTW_3D(temp,phi_k[indx1],FORWARD);
			}
		break;
	default:
		cout<<"Invalid elastic couple option!"<<endl;
		exit(1);
	}
	delete[] temp;

	return 0;
}

int Gelastic::Calc_Bpq()
{
	int g1,v1,g2,v2,i,j,k;
	long int indx11,indx2,indx1_1,indx1_2;
	int ii,jj,kk,ll;
	float det;
	float Bpp_ave;

	float omega[3][3],i_omega[3][3];

	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++)
	{
		indx2=index2(i,j,k);
		if(indx2!=0)
		{
//	Calc omega
			for(ii=0;ii<3;ii++)
			for(jj=0;jj<3;jj++)
			{
				i_omega[ii][jj]=0;
				for(kk=0;kk<3;kk++)
				for(ll=0;ll<3;ll++)
					i_omega[ii][jj]+=Cijkl[ii][kk][ll][jj]
							*(g[kk])[indx2]*(g[ll])[indx2]/g_sqr[indx2];
			}

			det=i_omega[0][0]*i_omega[1][1]*i_omega[2][2]
				+i_omega[0][1]*i_omega[1][2]*i_omega[2][0]
				+i_omega[1][0]*i_omega[2][1]*i_omega[0][2]
				-i_omega[2][0]*i_omega[0][2]*i_omega[1][1]
				-i_omega[0][1]*i_omega[1][0]*i_omega[2][2]
				-i_omega[2][1]*i_omega[1][2]*i_omega[0][0];

			omega[0][0]=(i_omega[1][1]*i_omega[2][2]-i_omega[1][2]*i_omega[2][1])/det;
			omega[0][1]=-(i_omega[1][0]*i_omega[2][2]-i_omega[1][2]*i_omega[2][0])/det;
			omega[0][2]=(i_omega[1][0]*i_omega[2][1]-i_omega[1][1]*i_omega[2][0])/det;
			omega[1][0]=-(i_omega[0][1]*i_omega[2][2]-i_omega[0][2]*i_omega[2][1])/det;
			omega[1][1]=(i_omega[0][0]*i_omega[2][2]-i_omega[0][2]*i_omega[2][0])/det;
			omega[1][2]=-(i_omega[0][0]*i_omega[2][1]-i_omega[0][1]*i_omega[2][0])/det;
			omega[2][0]=(i_omega[0][1]*i_omega[1][2]-i_omega[0][2]*i_omega[1][1])/det;
			omega[2][1]=-(i_omega[0][0]*i_omega[1][2]-i_omega[0][2]*i_omega[1][0])/det;
			omega[2][2]=(i_omega[0][0]*i_omega[1][1]-i_omega[0][1]*i_omega[1][0])/det;
//	Calc Bpq
			for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			for(g2=0;g2<ng;g2++)
			for(v2=0;v2<nv;v2++)
			{
				indx11=index11(g1,v1,g2,v2);
				indx1_1=index1(g1,v1);
				indx1_2=index1(g2,v2);

				(Bpq[indx11])[indx2]=0;

//				if(v1==10 || v2==10 ||v1==11||v2==11)	// Concentration Field
//					continue;

				for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				for(kk=0;kk<3;kk++)
				for(ll=0;ll<3;ll++)
				{
					(Bpq[indx11])[indx2]+=
								Cijkl[ii][jj][kk][ll]*(e0[indx1_1])[ii*3+jj]*(e0[indx1_2])[kk*3+ll]
								-(g[ii])[indx2]*(s0[indx1_1])[ii*3+jj]*omega[jj][kk]*(s0[indx1_2])[kk*3+ll]*(g[ll])[indx2]/g_sqr[indx2];
				}
			}
		}
		else // Bpq zero point
		{
			for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			for(g2=0;g2<ng;g2++)
			for(v2=0;v2<nv;v2++)
			{
				indx11=index11(g1,v1,g2,v2);
				indx1_1=index1(g1,v1);
				indx1_2=index1(g2,v2);

				switch(bc)
				{
				case FIX:
					(Bpq[indx11])[indx2]=0;

//					if(v1>3 || v2>3)	// Concentration Field
//						continue;

					for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					for(kk=0;kk<3;kk++)
					for(ll=0;ll<3;ll++)
					{
						(Bpq[indx11])[indx2]+=Cijkl[ii][jj][kk][ll]*(e0[indx1_1])[ii*3+jj]*(e0[indx1_2])[kk*3+ll];
					}
					break;
				case RELAX:
					(Bpq[indx11])[indx2]=0;
					break;
				default:
					cout<<"Invalid elastic boudary condition!"<<endl;
					exit(1);
				}
			}
		}
	}
   cd=0;
	float Bpq_min=10000,i_min=0,j_min=0,k_min=0;
	if(cd==1)
	{
		for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
	//				if(v1>3)	// Concentration Field
	//					continue;
			indx11=index11(g1,v1,g1,v1);
	//		indx11=index11(0,0,0,0);
			
			Bpp_ave=0;
			for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				Bpp_ave+=Bpq[indx11][indx2];
			}
			Bpp_ave/=nx*ny*nz;
			
			for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				Bpq[indx11][indx2]-=Bpp_ave;
			}
		}
//	cout<<endl<<"Average Bpq is:"<<Bpp_ave<<endl;
	}

//	Output Bpq for test			
	if(cd==2)
	{
		for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx11=index11(0,0,0,0);
			for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				if(Bpq[indx11][indx2]<=Bpq_min)
				{
					Bpq_min=Bpq[indx11][indx2];
				}
                       }
			for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				Bpq[indx11][indx2]-=Bpq_min;
			}
       }
       }
	
	ofstream fout("Bpq_minus.vtk",ios::out);
	Output_VTK_header(&fout,nx,ny,nz);
							
	indx11=index11(0,0,0,0);

	for(k=-nz/2;k<nz/2;k++) // ignore if for 2D
		for(j=-ny/2;j<ny/2;j++)
			for(i=-nx/2;i<nx/2;i++)
			{
				indx2=index2((i+nx)%nx,(j+ny)%ny,(k+nz)%nz);
				fout
//					<<i<<'\t'<<j<<'\t'<<k<<'\t'
					<<Bpq[indx11][indx2]<<endl;

				if(Bpq[indx11][indx2]<=Bpq_min && (i||j||k))
				{
					Bpq_min=Bpq[indx11][indx2];
					i_min=i;
					j_min=j;
					k_min=k;
				}
			}
	fout.flush();
	fout.close();

	cout<<endl<<"Minimux Bpq is:"<<endl<<i_min<<' '<<j_min<<' '<<k_min<<'\t'<<Bpq_min<<endl;


	return 0;
}

int Gelastic::Calc_Bpd()
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	int ii,jj,kk,ll,mm,nn;
	float det;

	float omega[3][3],i_omega[3][3];

	for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
	for(k=0;k<nz;k++)
	{
		indx2=index2(i,j,k);
		if(indx2!=0)
		{
//	Calc omega
			for(ii=0;ii<3;ii++)
			for(jj=0;jj<3;jj++)
			{
				i_omega[ii][jj]=0;
				for(kk=0;kk<3;kk++)
				for(ll=0;ll<3;ll++)
					i_omega[ii][jj]+=Cijkl[ii][kk][ll][jj]
									*(g[kk])[indx2]*(g[ll])[indx2]/g_sqr[indx2];
			}

			det=i_omega[0][0]*i_omega[1][1]*i_omega[2][2]
				+i_omega[0][1]*i_omega[1][2]*i_omega[2][0]
				+i_omega[1][0]*i_omega[2][1]*i_omega[0][2]
				-i_omega[2][0]*i_omega[0][2]*i_omega[1][1]
				-i_omega[0][1]*i_omega[1][0]*i_omega[2][2]
				-i_omega[2][1]*i_omega[1][2]*i_omega[0][0];

			omega[0][0]=(i_omega[1][1]*i_omega[2][2]-i_omega[1][2]*i_omega[2][1])/det;
			omega[0][1]=-(i_omega[1][0]*i_omega[2][2]-i_omega[1][2]*i_omega[2][0])/det;
			omega[0][2]=(i_omega[1][0]*i_omega[2][1]-i_omega[1][1]*i_omega[2][0])/det;
			omega[1][0]=-(i_omega[0][1]*i_omega[2][2]-i_omega[0][2]*i_omega[2][1])/det;
			omega[1][1]=(i_omega[0][0]*i_omega[2][2]-i_omega[0][2]*i_omega[2][0])/det;
			omega[1][2]=-(i_omega[0][0]*i_omega[2][1]-i_omega[0][1]*i_omega[2][0])/det;
			omega[2][0]=(i_omega[0][1]*i_omega[1][2]-i_omega[0][2]*i_omega[1][1])/det;
			omega[2][1]=-(i_omega[0][0]*i_omega[1][2]-i_omega[0][2]*i_omega[1][0])/det;
			omega[2][2]=(i_omega[0][0]*i_omega[1][1]-i_omega[0][1]*i_omega[1][0])/det;
//	Calc Bpd
			for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);

				for(mm=0;mm<3;mm++)
				for(nn=0;nn<3;nn++)
					((Bpd[mm][nn])[indx1])[indx2]=0;

//				if(v1>3)	// Concentration Field
//					continue;

				for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				for(kk=0;kk<3;kk++)
				for(ll=0;ll<3;ll++)
				{
					for(mm=0;mm<3;mm++)
					for(nn=0;nn<3;nn++)
					{
						((Bpd[mm][nn])[indx1])[indx2]+=
									Cijkl[ii][jj][mm][nn]*(e0[indx1])[ii*3+jj]/9.0
									-(g[ii])[indx2]*(s0[indx1])[ii*3+jj]*omega[jj][kk]*Cijkl[kk][ll][mm][nn]*(g[ll])[indx2]/g_sqr[indx2]; // first term sum of ii,jj
					}
				}
			}
		}
		else // Bpd zero point
		{
			for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);

				switch(bc)
				{
				case FIX:
					for(mm=0;mm<3;mm++)
					for(nn=0;nn<3;nn++)
						((Bpd[mm][nn])[indx1])[indx2]=0;

//					if(v1>3)	// Concentration Field
//						continue;

					for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					for(kk=0;kk<3;kk++)
					for(ll=0;ll<3;ll++)
					{
						for(mm=0;mm<3;mm++)
						for(nn=0;nn<3;nn++)
							((Bpd[mm][nn])[indx1])[indx2]+=Cijkl[ii][jj][mm][nn]*(e0[indx1])[ii*3+jj]/9.0; // summation of ii,jj
					}
					break;
				case RELAX:
					for(mm=0;mm<3;mm++)
					for(nn=0;nn<3;nn++)
						((Bpd[mm][nn])[indx1])[indx2]=0;
					break;
				case CONFIG_DEPT:
					cout<<"Not availiable boudary condition!"<<endl;
					exit(1);
					break;
				default:
					cout<<"Invalid elastic boudary condition!"<<endl;
					exit(1);
				}
			}
		}
	}
		if(cd==1)
			cout<<"Config-dept option is invalid for Bpd calculation!"<<endl;			;
	return 0;
}

int Gelastic::Set_g()
{
	int i,j,k;
	long int indx2;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				if(i<nx/2.0)
					(g[0])[indx2]=2*Pi/nx*i;
				else
					(g[0])[indx2]=2*Pi/nx*(i-nx);

				if(j<ny/2.0)
					(g[1])[indx2]=2*Pi/ny*j;
				else
					(g[1])[indx2]=2*Pi/ny*(j-ny);
				
				if(k<nz/2.0)
					(g[2])[indx2]=2*Pi/nz*k;
				else
					(g[2])[indx2]=2*Pi/nz*(k-nz);

				g_sqr[indx2]=(g[0])[indx2]*(g[0])[indx2]
						+(g[1])[indx2]*(g[1])[indx2]
						+(g[2])[indx2]*(g[2])[indx2];
			}

	return 0;
}

int Gelastic::AppliedStress(float s_applied[][3])
{
	int g1,v1;
	long int indx1,indx2;
	int mm,nn,i,j,k,mnindx;
	ofstream *fout,vout,eout;

	if(bc!=RELAX)
	{
		cout<<"Invalid boundary condition for applied stress!"<<endl;
		exit(1);
	}

	float **stress=new float*[3*3];
	assert(stress!=NULL);
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			mnindx=mm*3+nn;
			stress[mnindx]=new float[nx*ny*nz];
			assert(stress[mnindx]!=NULL);
		}
	//Load_StressField(stress);
        float **potential= new float*[ng*nv];
	fout=new ofstream[ng*nv];
	assert(fout!=NULL);
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			potential[indx1]=new float[nx*ny*nz];
			assert(potential[indx1]!=NULL);
		}
	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						potential[indx1][indx2]=0;
						for(mm=0;mm<3;mm++)
							for(nn=0;nn<3;nn++)
							{
								mnindx=mm*3+nn;
								potential[indx1][indx2]=stress[mnindx][indx2]*e0[indx1][mnindx];
							}
							// InterEnergy for each variant
					//		fout[indx1]<<-energy[indx1]<<endl;
					}
			}


	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
			s_app[mm][nn]=s_applied[mm][nn]*s_app_mag;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			s_app_miu[indx1]=0;
			for(mm=0;mm<3;mm++)
				for(nn=0;nn<3;nn++)
					s_app_miu[indx1]-=s_app[mm][nn]*(e0[indx1])[mm*3+nn];
	       cout<<"app_miu="<<s_app_miu[indx1]<<endl;
		}
	return 0;
}
int Gelastic::UpdateStress(float s_applied[][3])
{
	int g1,v1;
	long int indx1;
	int mm,nn;
	int num=0;
       // num++;
	if(bc!=RELAX)
	{
		cout<<"Invalid boundary condition for applied stress!"<<endl;
		exit(1);
	}
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
			s_app[mm][nn]+=s_applied[mm][nn]*(s_app_inc);
		//	s_app[mm][nn]+=s_applied[mm][nn]*(-10);
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			s_app_miu[indx1]=0;
			for(mm=0;mm<3;mm++)
				for(nn=0;nn<3;nn++)
					s_app_miu[indx1]-=s_app[mm][nn]*(e0[indx1])[mm*3+nn];
	       cout<<"update_app_miu="<<s_app_miu[indx1]<<endl;
		}
	return 0;


}
int Gelastic::AppliedStrain(float e_applied[][3])
{
	int g1,v1;
	long int indx1;
	int ii,jj,kk,ll;

	if(bc!=FIX)
	{
		cout<<"Invalid boundary condition for applied strain!"<<endl;
		exit(1);
	}
	for(ii=0;ii<3;ii++)
		for(jj=0;jj<3;jj++)
			e_app[ii][jj]=e_applied[ii][jj]*e_app_mag;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			e_app_miu[indx1]=0;
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
					for(kk=0;kk<3;kk++)
						for(ll=0;ll<3;ll++)
							//e_app_miu[indx1]-=1.1*Cijkl[ii][jj][kk][ll]*e_app[ii][jj]*e0[indx1][kk*3+ll];
							e_app_miu[indx1]-=Cijkl[ii][jj][kk][ll]*e_app[ii][jj]*e0[indx1][kk*3+ll];
	       cout<<"app_miu="<<e_app_miu[indx1]<<endl;
		}
	return 0;
}
int Gelastic::UpdateStrain(float e_applied[][3])
{
	int g1,v1;
	long int indx1;
	int mm,nn,ii,jj,kk,ll;

	if(bc!=FIX)
	{
		cout<<"Invalid boundary condition for applied stress!"<<endl;
		exit(1);
	}
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
			e_app[mm][nn]=e_applied[mm][nn]*(e_app_inc);
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
	//		s_app_miu[indx1]=0;
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
					for(kk=0;kk<3;kk++)
						for(ll=0;ll<3;ll++)
							e_app_miu[indx1]-=Cijkl[ii][jj][kk][ll]*e_app[ii][jj]*(e0[indx1])[kk*3+ll];
	       cout<<"app_miu="<<e_app_miu[indx1]<<endl;
		}
	return 0;


}

int Gelastic::Output_Stress(float **eta,int *gs) // Total stress field: including applied stress/strain
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	int mm,nn,kk,ll;
	char f0[20], f1[20], f2[20], f3[20];
	char f4[20], f5[20], f6[20], f7[20], f8[20];

//	if(output_stress!=1 || Bpd==NULL)
	if(Bpd==NULL)
	{
		cout<<"Invalid input for output stress!"<<endl;
		exit(1);
	}

	Calc_phi_k(eta);

	fftw_complex *miu_k;
	miu_k=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*ny*nz);
	assert(miu_k!=NULL);
	float *miu=new float[nx*ny*nz];
	assert(miu!=NULL);

	char fn[]="stress0.vtk";

	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);

						miu_k[indx2][0]=0;
						miu_k[indx2][1]=0;
						
						for(g1=0;g1<ng;g1++)
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);

							miu_k[indx2][0]-=((Bpd[mm][nn])[indx1])[indx2]*(phi_k[indx1])[indx2][0];
							miu_k[indx2][1]-=((Bpd[mm][nn])[indx1])[indx2]*(phi_k[indx1])[indx2][1];
						}				
					}
			gfft.FFTW_3D(miu_k,miu,BACKWARD);
			
		//	fn[6]=(char)((int)'0'+mm*3+nn);
		//	ofstream fout(fn,ios::out);
		//	Output_VTK_header(&fout,nx,ny,nz);
		//	fout<<fixed;
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						
						miu[indx2]*=scale;

						if(fabs(s_app_mag)>1E-5)
							miu[indx2]+=s_app[mm][nn];

						if(fabs(e_app_mag)>1E-5)
						{
							for(kk=0;kk<3;kk++)
								for(ll=0;ll<3;ll++)
								{
									miu[indx2]+=Cijkl[mm][nn][kk][ll]*e_app[kk][ll];
								}
						}

		//				fout<<setprecision(9)<<miu[indx2]<<endl;
	      				s_field[mm*3+nn][indx2]=miu[indx2]/100;
	      			}
//			fout.flush();
//			fout.close();
		}
	
	fftw_free(miu_k);
	delete[] miu;

	return 0;
}


int Gelastic::Output_Strain(float **eta,int *gs, int out) // Total stress field: including applied stress/strain
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	int mm,nn,ii,jj,kk,ll;
//FILE *fp=NULL;
	char f0[20], f1[20], f2[20], f3[20];
	char f4[20], f5[20], f6[20], f7[20], f8[20];
	static int num=0;
	float e0_total[3][3]={0},s_total[3][3]={0};
	float e1,fraction;

/*	if(output_stress!=1 || Bpd==NULL)
	{
		cout<<"Invalid input for output stress!"<<endl;
		exit(1);
	}
*/
	float **e_field=new float *[9];
	assert(e_field!=NULL);
	for(int i=0;i<9;i++)
	{
		e_field[i]=new float [nx*ny*nz];
		assert(e_field[i]!=NULL);
	}

	Set_Sijkl_Hphase();
/*
				for(mm=0;mm<3;mm++)
					for(nn=0;nn<3;nn++)
						for(kk=0;kk<3;kk++)
							for(ll=0;ll<3;ll++)
							{
	printf("S_%d_%d_%d_%d=%lf\n",mm,nn,kk,ll,Sijkl[mm][nn][kk][ll]);
               }
*/	   
/*
        sprintf(f0,"strain11_%02d.vtk",num);
	sprintf(f1,"strain12_%02d.vtk",num);
	sprintf(f2,"strain13_%02d.vtk",num);
	sprintf(f3,"strain21_%02d.vtk",num);
        sprintf(f4,"strain22_%02d.vtk",num);
	sprintf(f5,"strain23_%02d.vtk",num);
	sprintf(f6,"strain31_%02d.vtk",num);
	sprintf(f7,"strain32_%02d.vtk",num);
	sprintf(f8,"strain33_%02d.vtk",num);

        num++;
	ofstream f0_out(f0,ios::out);
	ofstream f1_out(f1,ios::out);
	ofstream f2_out(f2,ios::out);
	ofstream f3_out(f3,ios::out);
	ofstream f4_out(f4,ios::out);
	ofstream f5_out(f5,ios::out);
	ofstream f6_out(f6,ios::out);
	ofstream f7_out(f7,ios::out);
	ofstream f8_out(f8,ios::out);
	Output_VTK_header(&f0_out,nx,ny,nz);
	Output_VTK_header(&f1_out,nx,ny,nz);
	Output_VTK_header(&f2_out,nx,ny,nz);
	Output_VTK_header(&f3_out,nx,ny,nz);
	Output_VTK_header(&f4_out,nx,ny,nz);
	Output_VTK_header(&f5_out,nx,ny,nz);
	Output_VTK_header(&f6_out,nx,ny,nz);
	Output_VTK_header(&f7_out,nx,ny,nz);
	Output_VTK_header(&f8_out,nx,ny,nz);
*/

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(mm=0;mm<3;mm++)
					for(nn=0;nn<3;nn++)
					{
				//		e_field[mm*3+nn][indx2]=0;
						for(kk=0;kk<3;kk++)
							for(ll=0;ll<3;ll++)
							{
								e_field[mm*3+nn][indx2]+=Sijkl[mm][nn][kk][ll]*s_field[kk*3+ll][indx2];
							}
						float tp=0;
						for(int v1=0;v1<nv;v1++)
						{
							tp+=eta[v1][indx2]*e0[v1][mm*3+nn];
						}

						e_field[mm*3+nn][indx2]+=tp;
					}
			}


//	for(mm=0;mm<3;mm++)
//		for(nn=0;nn<3;nn++)
//		{
//			fn[6]=(char)((int)'0'+mm*3+nn);
		//	ofstream fout(fn,ios::out);
		//	Output_VTK_header(&fout,nx,ny,nz);
//			fout<<fixed;
/*			
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						
                 	for(mm=0;mm<3;mm++)
		               for(nn=0;nn<3;nn++)
	 	                     { 
						float tp=0;
						for(int v1=0;v1<nv;v1++)
						{
							tp+=eta[v1][indx2]*e0[v1][mm*3+nn];
						}

						e_field[mm*3+nn][indx2]+=tp;
						}
				//	fprintf(fp,"%.6lf \n",e_field[mm*3+nn][indx2]);	
//	fprintf(fp,"%.6lf \n",e_field[mm*3+nn][indx2]);	
*/
/*						f0_out<<setprecision(9)<<e_field[0][indx2]<<endl;
						f1_out<<setprecision(9)<<e_field[1][indx2]<<endl;
						f2_out<<setprecision(9)<<e_field[2][indx2]<<endl;
						f3_out<<setprecision(9)<<e_field[3][indx2]<<endl;
						f4_out<<setprecision(9)<<e_field[4][indx2]<<endl;
						f5_out<<setprecision(9)<<e_field[5][indx2]<<endl;
						f6_out<<setprecision(9)<<e_field[6][indx2]<<endl;
						f7_out<<setprecision(9)<<e_field[7][indx2]<<endl;
						f8_out<<setprecision(9)<<e_field[8][indx2]<<endl;
*/				
//				}//end of space loop
	//		fout.flush();
	//		fout.close()
	//		}
/*	f0_out.flush();
	f0_out.close();
	f1_out.flush();
	f1_out.close();
	f2_out.flush();
	f2_out.close();
	f3_out.flush();
	f3_out.close();
	f4_out.flush();
	f4_out.close();

	f5_out.flush();
	f5_out.close();
	
	f6_out.flush();
	f6_out.close();
	f7_out.flush();
	f7_out.close();
	f8_out.flush();
	f8_out.close();
*/

//	ofstream fout("Total Strain.dat",ios::out);
//	cout<<"Total Strain: "<<num<<endl;

	cout<<"Calculate average total strain!"<<endl;
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			e0_total[mm][nn]=0;

			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);
						for(g1=0;g1<ng;g1++)
						{
							if(g1!=gs[indx2])
								continue;

							for(v1=0;v1<nv;v1++)
							{
								indx1=index1(g1,v1);
								e1=(eta[indx1])[indx2];

								e0_total[mm][nn]+=e_field[mm*3+nn][indx2]/nx/ny/nz;
							}
						}
					}
		//	cout<<e0_total[mm][nn]<<endl;
		}
			cout<<"Average Strain11="<<e0_total[0][0]<<endl;
//	fout.flush();
//	fout.close();

	
	return 0;
}

int Gelastic::Set_Sijkl_Hphase()
{
	int ii,jj,kk,ll;
        voigt66 Cij0,Sij0;

	//	Sijkl for rotated NiTi B2
	cout<<"Strain field is only for rotation NiTi B2!"<<endl;
//	chg_basis_Kelvin_4(Cijkl,Cij0);
      C4_loop{
            int p = (mi+1)*(mi==mj)+(1-(mi==mj))*(7-mi-mj);
            int q = (mk+1)*(mk==ml)+(1-(mk==ml))*(7-mk-ml);
            real t1 = ((real)(mi==mj)+(1.-(mi==mj))/RSQ2);
            real t2 = ((real)(mk==ml)+(1.-(mk==ml))/RSQ2);
            Cij0[p-1][q-1] = t1*t2*Cijkl[mi][mj][mk][ml];
        }
	C6_loop{
        Sij0[mi][mj] = Cij0[mi][mj];
    	}
    	LU_inv_66(Sij0);
//	chg_basis_Kelvin_3(Sij0,Sijkl);
        C4_loop{
            int p = (mi+1)*(mi==mj)+(1-(mi==mj))*(7-mi-mj);
            int q = (mk+1)*(mk==ml)+(1-(mk==ml))*(7-mk-ml);
            real t1 = ((real)(mi==mj)+(1.-(mi==mj))/RSQ2);
            real t2 = ((real)(mk==ml)+(1.-(mk==ml))/RSQ2);
            Sijkl[mi][mj][mk][ml] = Sij0[p-1][q-1]/t1/t2;
        }

        C4_loop{
	
            Sijkl[mi][mj][mk][ml] = Sijkl[mi][mj][mk][ml]/scale;

	}
/*
	for(ii=0;ii<3;ii++)
		for(jj=0;jj<3;jj++)
			for(kk=0;kk<3;kk++)
				for(ll=0;ll<3;ll++)
				{
					Sijkl[ii][jj][kk][ll]=0;
				}
	Sijkl[0][0][0][0]=1.8720;
	Sijkl[1][1][1][1]=1.8720;
	Sijkl[2][2][2][2]=1.8720;
	
	Sijkl[0][0][1][1]=Sijkl[1][1][0][0]=-0.8307;
	Sijkl[1][1][2][2]=Sijkl[2][2][1][1]=-0.8307;
	Sijkl[2][2][0][0]=Sijkl[0][0][2][2]=-0.8307;

	Sijkl[0][1][0][1]=Sijkl[0][1][1][0]=Sijkl[1][0][0][1]=Sijkl[1][0][1][0]=2.1739;
	Sijkl[1][2][1][2]=Sijkl[1][2][2][1]=Sijkl[2][1][1][2]=Sijkl[2][1][2][1]=2.1739;
	Sijkl[2][0][2][0]=Sijkl[2][0][0][2]=Sijkl[0][2][2][0]=Sijkl[0][2][0][2]=5.4054;

	Sijkl[0][1][0][1]=Sijkl[0][1][1][0]=Sijkl[1][0][0][1]=Sijkl[1][0][1][0]=0.5435;
	Sijkl[1][2][1][2]=Sijkl[1][2][2][1]=Sijkl[2][1][1][2]=Sijkl[2][1][2][1]=0.5435;
	Sijkl[2][0][2][0]=Sijkl[2][0][0][2]=Sijkl[0][2][2][0]=Sijkl[0][2][0][2]=0.5435;
*/
return 0;
}


int Gelastic::Output_InterEnergy() // -stress*strain
{
	int g1,v1,i,j,k;
	long indx1,indx2,mnindx;
	int mm,nn;
	float *energy,max_energy=0;
	int maxID=-1;
	ofstream *fout,vout,eout;

	float highest_energy=0;
	ofstream hout;
	
	if(output_interenergy!=1)
	{
		cout<<"Invalid input for output interaction fields!"<<endl;
		exit(1);
	}

//	Initialize stress field
	float **stress=new float*[3*3];
	assert(stress!=NULL);
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			mnindx=mm*3+nn;
			stress[mnindx]=new float[nx*ny*nz];
			assert(stress[mnindx]!=NULL);
		}
//	Load_StressField(stress);

//	Calculate Interaction Energy
	energy=new float[ng*nv];
	assert(energy!=NULL);
	fout=new ofstream[ng*nv];
	assert(fout!=NULL);

	char fn1[]="InterEnergy0.vtk";
	char fn2[]="MaxEnergyID.vtk";
	char fn3[]="MaxEnergy.vtk";

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			fn1[11]=(char)((int)'0'+indx1);
			fout[indx1].open(fn1,ios::out);
			Output_VTK_header(&fout[indx1],nx,ny,nz);
		}
	vout.open(fn2,ios::out);
	Output_VTK_header(&vout,nx,ny,nz);
	eout.open(fn3,ios::out);
	Output_VTK_header(&eout,nx,ny,nz);

	hout.open("LowestInterEnergy.txt",ios::out); // Record Lowest Interaction Energy

	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						energy[indx1]=0;
						for(mm=0;mm<3;mm++)
							for(nn=0;nn<3;nn++)
							{
								mnindx=mm*3+nn;
								energy[indx1]+=stress[mnindx][indx2]*e0[indx1][mnindx];
							}
							// InterEnergy for each variant
							fout[indx1]<<-energy[indx1]<<endl;
					}
	// max InterEnergy variant ID
				maxID=0;
				max_energy=energy[0];
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						if(energy[indx1]>max_energy)
						{
							max_energy=energy[indx1];
							maxID=indx1;			
						}
					}
				// max InterEnergy ID
				if(max_energy>0.00001)
					vout<<maxID<<endl;
				else 
					vout<<-1<<endl;
				// max InterEnergy value (negative value)
				eout<<-max_energy<<endl;

				highest_energy=highest_energy>max_energy?highest_energy:max_energy;
			}

	hout<<"Lowest Interaction Energy = "<<-highest_energy<<endl;
	hout<<-highest_energy*2.84<<" J/mm3"<<endl;
	hout<<-highest_energy*2.84*1.771*10/2.0<<" KJ/mol atom"<<endl;

	vout.flush();
	vout.close();
	eout.flush();
	eout.close();
	hout.flush();
	hout.close();

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			fout[indx1].flush();
			fout[indx1].close();
		}

	delete[] fout;
	delete[] energy;
//	Release stress field
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			mnindx=mm*3+nn;
			delete[] stress[mnindx];
			stress[mnindx]=NULL;
		}
	delete[] stress;
	stress=NULL;
	
	
	return 0;
}

int Gelastic::Load_StressField(float **stress)
{
	
	int i,j,k;
	long int indx2,mnindx;
	int mm,nn;
	char temp[80]={0};
	char fn[]="stress0.vtk";
	ifstream fin;
        cout<<"Stress is loaded!"<<endl;
	if(stress==NULL)
	{
		cout<<"No space for stress loading!"<<endl;
		exit(1);
	}
	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			mnindx=mm*3+nn; // component index
			fn[6]=(char)((int)'0'+mnindx);
			fin.open(fn,ios::in);
			if(!fin)
			{
				cerr<<"stress input can not be found!"<<endl;
				exit(1);
			}
			for(i=0;i<30;i++) // ignore VTK header
				fin>>temp;
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						fin>>stress[mnindx][indx2];
					}
			fin.close();
		}

	// s11+s22+s33
	ofstream fout("stress-sum.vtk",ios::out);
	Output_VTK_header(&fout,nx,ny,nz);
	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
			{
				indx2=index2(i,j,k);
				fout<<stress[0][indx2]+stress[4][indx2]+stress[8][indx2]<<endl;
			}
	fout.flush();
	fout.close();
	
	return 0;
}

int Gelastic::Reset_SFTS_MT() // for InteractionEnergy calculation in NiTiPt
{
	int g1,v1;
	int ii,jj,kk;
	long int indx1,iindx2;
	float aa[3]={1.0192,0.8929,1.0789};

	float eye[3][3]={0};
	eye[0][0]=eye[1][1]=eye[2][2]=1;


//	Set U
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			switch(v1)
			{
			case 0:
				// Orthorhombic
				(u[indx1])[0]=(aa[0]+aa[2])/2.0;
				(u[indx1])[1]=0;
				(u[indx1])[2]=(aa[0]-aa[2])/2.0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=aa[1];
				(u[indx1])[5]=0;
				(u[indx1])[6]=(aa[0]-aa[2])/2.0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=(aa[0]+aa[2])/2.0;

				break;

			case 1:
				// Orthorhombic
				(u[indx1])[0]=(aa[0]+aa[2])/2.0;
				(u[indx1])[1]=0;
				(u[indx1])[2]=-(aa[0]-aa[2])/2.0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=aa[1];
				(u[indx1])[5]=0;
				(u[indx1])[6]=-(aa[0]-aa[2])/2.0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=(aa[0]+aa[2])/2.0;

				break;

			case 2:
				// Orthorhombic
				(u[indx1])[0]=(aa[0]+aa[2])/2.0;
				(u[indx1])[1]=(aa[0]-aa[2])/2.0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=(aa[0]-aa[2])/2.0;
				(u[indx1])[4]=(aa[0]+aa[2])/2.0;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=aa[1];

				break;

			case 3:
				// Orthorhombic
				(u[indx1])[0]=(aa[0]+aa[2])/2.0;
				(u[indx1])[1]=-(aa[0]-aa[2])/2.0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=-(aa[0]-aa[2])/2.0;
				(u[indx1])[4]=(aa[0]+aa[2])/2.0;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=aa[1];

				break;

			case 4:
				// Orthorhombic
				(u[indx1])[0]=aa[1];
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=(aa[0]+aa[2])/2.0;
				(u[indx1])[5]=(aa[0]-aa[2])/2.0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=(aa[0]-aa[2])/2.0;
				(u[indx1])[8]=(aa[0]+aa[2])/2.0;

				break;

			case 5:
				// Orthorhombic
				(u[indx1])[0]=aa[1];
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=(aa[0]+aa[2])/2.0;
				(u[indx1])[5]=-(aa[0]-aa[2])/2.0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=-(aa[0]-aa[2])/2.0;
				(u[indx1])[8]=(aa[0]+aa[2])/2.0;

				break;

			default:
				(u[indx1])[0]=1;
				(u[indx1])[1]=0;
				(u[indx1])[2]=0;
				(u[indx1])[3]=0;
				(u[indx1])[4]=1;
				(u[indx1])[5]=0;
				(u[indx1])[6]=0;
				(u[indx1])[7]=0;
				(u[indx1])[8]=1;				
				break;
			}
		}

// Set SFTS
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
				for(jj=0;jj<3;jj++)
				{
					iindx2=ii*3+jj;
					(e0[indx1])[iindx2]=0;
					for(kk=0;kk<3;kk++)
						e0[indx1][iindx2]+=(u[indx1])[kk*3+ii]*(u[indx1])[kk*3+jj];

					e0[indx1][iindx2]-=eye[ii][jj];
					e0[indx1][iindx2]/=2.0;
				}
		}

// Output SFTS for InterEnergy Calculation

	cout<<"Rest SFTS for InterEnergy Calculation:"<<endl;
	for(g1=0;g1<ng;g1++)
	{
		cout<<"In grain "<<g1<<endl;
		for(v1=0;v1<nv;v1++)
		{
			cout<<"variant "<<v1<<endl;
			indx1=index1(g1,v1);
			for(ii=0;ii<3;ii++)
			{
				for(jj=0;jj<3;jj++)
					cout<<e0[indx1][ii*3+jj]<<'\t';
				cout<<endl;
			}
			cout<<endl;
		}
		cout<<endl;
	}


	return 0;
}

int Gelastic::Output_TranStrain(float **eta,int *gs)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	int mm,nn,kk,ll;
	float e1,fraction;
	float e0_total[3][3]={0},s_total[3][3]={0};

	ofstream fout("TranStrain.dat",ios::out);
	fout<<"Total TranStrain: "<<endl;

	for(mm=0;mm<3;mm++)
		for(nn=0;nn<3;nn++)
		{
			e0_total[mm][nn]=0;

			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);
						for(g1=0;g1<ng;g1++)
						{
							if(g1!=gs[indx2])
								continue;

							for(v1=0;v1<nv;v1++)
							{
								indx1=index1(g1,v1);
								e1=(eta[indx1])[indx2];

								e0_total[mm][nn]+=e1*(e0[indx1])[mm*3+nn]/nx/ny/nz;
							}
						}
					}
			fout<<e0_total[mm][nn]<<endl;
		}

	if(bc==FIX)
	{
		fout<<endl<<"Total Stress for strain controlled condition: "<<endl;
		for(mm=0;mm<3;mm++)
			for(nn=0;nn<3;nn++)
			{
				s_total[mm][nn]=0;
				for(kk=0;kk<3;kk++)
					for(ll=0;ll<3;ll++)
					{
						s_total[mm][nn]+=Cijkl[mm][nn][kk][ll]*(e_app[kk][ll]-e0_total[kk][ll]);
					}
				fout<<s_total[mm][nn]<<endl;
			}
	}

	fout<<endl<<"Fraction of each variant: "<<endl;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			fraction=0;
			indx1=index1(g1,v1);
			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);
						fraction+=(eta[indx1])[indx2];
					}
			fout<<fraction/nx/ny/nz<<endl;
		}	

	fout.flush();
	fout.close();

	return 0;
}



int Gelastic::Output_VTK_header(ofstream *p_fout,int nxx,int nyy,int nzz)
{
	*p_fout<<"# vtk DataFile Version 2.0"<<endl;
	*p_fout<<"Volume example"<<endl<<"ASCII"<<endl<<"DATASET STRUCTURED_POINTS"<<endl;
	//*p_fout<<"DIMENSIONS"<<'\t'<<nxx<<'\t'<<nyy<<'\t'<<nzz<<endl;
	*p_fout<<"DIMENSIONS"<<" "<<nxx<<" "<<nyy<<" "<<nzz<<endl;
	*p_fout<<"ASPECT_RATIO 1 1 1"<<endl;
	*p_fout<<"ORIGIN 0 0 0"<<endl;
	*p_fout<<"POINT_DATA"<<" "<<nxx*nyy*nzz<<endl;
	*p_fout<<"SCALARS volume_scalars float 1"<<endl;
	*p_fout<<"LOOKUP_TABLE default"<<endl;

	return 0;
}

int Gelastic::Set_SFTS_NiPtTi()  // Used for both chemical and structural dependence of SFTS
{
	int g1,v1,ii,jj;
	long int indx1;

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			switch(v1)
			{
			case 0:
				e0[indx1][0]=-0.0126;
				e0[indx1][1]=-0.0016;
				e0[indx1][2]=-0.0016;
				e0[indx1][3]=-0.0016;
				e0[indx1][4]=-0.0126;
				e0[indx1][5]=-0.0016;
				e0[indx1][6]=-0.0016;
				e0[indx1][7]=-0.0016;
				e0[indx1][8]=-0.0126;

				break;

			case 1:
				e0[indx1][0]=-0.0126;
				e0[indx1][1]=0.0016;
				e0[indx1][2]=0.0016;
				e0[indx1][3]=0.0016;
				e0[indx1][4]=-0.0126;
				e0[indx1][5]=-0.0016;
				e0[indx1][6]=0.0016;
				e0[indx1][7]=-0.0016;
				e0[indx1][8]=-0.0126;

				break;

			case 2:
				e0[indx1][0]=-0.0126;
				e0[indx1][1]=0.0016;
				e0[indx1][2]=-0.0016;
				e0[indx1][3]=0.0016;
				e0[indx1][4]=-0.0126;
				e0[indx1][5]=0.0016;
				e0[indx1][6]=-0.0016;
				e0[indx1][7]=0.0016;
				e0[indx1][8]=-0.0126;

				break;

			case 3:
				e0[indx1][0]=-0.0126;
				e0[indx1][1]=-0.0016;
				e0[indx1][2]=0.0016;
				e0[indx1][3]=-0.0016;
				e0[indx1][4]=-0.0126;
				e0[indx1][5]=0.0016;
				e0[indx1][6]=0.0016;
				e0[indx1][7]=0.0016;
				e0[indx1][8]=-0.0126;

				break;

			case 4:
				for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					{
						if(ii!=jj)
							e0[indx1][ii*3+jj]=0;
						else
							e0[indx1][ii*3+jj]=-0.092;
					}
				break;
			case 5:
				for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					{
						if(ii!=jj)
							e0[indx1][ii*3+jj]=0;
						else
							e0[indx1][ii*3+jj]=-0.013;
					}
				break;
			default:
				break;
			}
		}
	return 0;
}

int Gelastic::Set_SFTS_PLAMINATE() // Laminate structure for P phase // Only variant 0/1/2
{
	int g1,v1,ii,jj;
	long int indx1;

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			switch(v1)
			{
			case 0:
				e0[indx1][0]=0.0003;
				e0[indx1][1]=-0.0005;
				e0[indx1][2]=-0.0022;
				e0[indx1][3]=-0.0005;
				e0[indx1][4]=0.0003;
				e0[indx1][5]=-0.0022;
				e0[indx1][6]=-0.0022;
				e0[indx1][7]=-0.0022;
				e0[indx1][8]=-0.0050;

				break;

			case 1:
				e0[indx1][0]=-0.0050;
				e0[indx1][1]=-0.0022;
				e0[indx1][2]=-0.0022;
				e0[indx1][3]=-0.0022;
				e0[indx1][4]=0.0003;
				e0[indx1][5]=-0.0005;
				e0[indx1][6]=-0.0022;
				e0[indx1][7]=-0.0005;
				e0[indx1][8]=0.0003;

				break;

			case 2:
				e0[indx1][0]=0.0003;
				e0[indx1][1]=-0.0022;
				e0[indx1][2]=-0.0005;
				e0[indx1][3]=-0.0022;
				e0[indx1][4]=-0.0050;
				e0[indx1][5]=-0.0022;
				e0[indx1][6]=-0.0005;
				e0[indx1][7]=-0.0022;
				e0[indx1][8]=0.0003;

				break;

			case 3:
				e0[indx1][0]=-0.0015;
				e0[indx1][1]=-0.0017;
				e0[indx1][2]=0.0017;
				e0[indx1][3]=-0.0017;
				e0[indx1][4]=-0.0015;
				e0[indx1][5]=0.0017;
				e0[indx1][6]=0.0017;
				e0[indx1][7]=0.0017;
				e0[indx1][8]=-0.0015;

				break;

			case 4:
				for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					{
						if(ii!=jj)
							e0[indx1][ii*3+jj]=0;
						else
							e0[indx1][ii*3+jj]=-0.1509;
					}
				break;
			case 5:
				for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					{
						if(ii!=jj)
							e0[indx1][ii*3+jj]=0;
						else
							e0[indx1][ii*3+jj]=0.089;
					}
				break;
			default:
				break;
			}
		}
	return 0;
}

int Gelastic::Set_SFTS_TiAlV()  // Used for TiAl (test for Rongpei)
{
	e0[0][0]=-0.051;
	e0[0][1]=0.023;
	e0[0][2]=0;
	e0[0][3]=0.023;
	e0[0][4]=0.034;
	e0[0][5]=0;
	e0[0][6]=0;
	e0[0][7]=0;
	e0[0][8]=-0.00029;

	return 0;
}

int Gelastic::Set_SFTS_NiTiHf()  // Used for H-phase
{
	int g1,v1,ii,jj;
	long int indx1;

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			switch(v1)
			{
			case 0:
				e0[indx1][0]=0.02291;
				e0[indx1][1]=0;
				e0[indx1][2]=0;
				e0[indx1][3]=0;
				e0[indx1][4]=0.00195;
				e0[indx1][5]=0.00726;
				e0[indx1][6]=0;
				e0[indx1][7]=0.00726;
				e0[indx1][8]=0.00195;

				break;

			case 1:
				e0[indx1][0]=0.02291;
				e0[indx1][1]=0;
				e0[indx1][2]=0;
				e0[indx1][3]=0;
				e0[indx1][4]=0.00195;
				e0[indx1][5]=-0.00726;
				e0[indx1][6]=0;
				e0[indx1][7]=-0.00726;
				e0[indx1][8]=0.00195;

				break;

			case 2:
				e0[indx1][0]=0.00195;
				e0[indx1][1]=0;
				e0[indx1][2]=0.00726;
				e0[indx1][3]=0;
				e0[indx1][4]=0.02291;
				e0[indx1][5]=0;
				e0[indx1][6]=0.00726;
				e0[indx1][7]=0;
				e0[indx1][8]=0.00195;

				break;

			case 3:
				e0[indx1][0]=0.00195;
				e0[indx1][1]=0;
				e0[indx1][2]=-0.00726;
				e0[indx1][3]=0;
				e0[indx1][4]=0.02291;
				e0[indx1][5]=0;
				e0[indx1][6]=-0.00726;
				e0[indx1][7]=0;
				e0[indx1][8]=0.00195;

				break;

			case 4:
				e0[indx1][0]=0.00195;
				e0[indx1][1]=0.00726;
				e0[indx1][2]=0;
				e0[indx1][3]=0.00726;
				e0[indx1][4]=0.00195;
				e0[indx1][5]=0;
				e0[indx1][6]=0;
				e0[indx1][7]=0;
				e0[indx1][8]=0.02291;

				break;
			case 5:
				e0[indx1][0]=0.00195;
				e0[indx1][1]=-0.00726;
				e0[indx1][2]=0;
				e0[indx1][3]=-0.00726;
				e0[indx1][4]=0.00195;
				e0[indx1][5]=0;
				e0[indx1][6]=0;
				e0[indx1][7]=0;
				e0[indx1][8]=0.02291;

				break;
				
			case 6:
				for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					{
						if(ii!=jj)
							e0[indx1][ii*3+jj]=0;
						else
							e0[indx1][ii*3+jj]=0;
					}
				break;

			case 7:
				for(ii=0;ii<3;ii++)
					for(jj=0;jj<3;jj++)
					{
						if(ii!=jj)
							e0[indx1][ii*3+jj]=0;
						else
							e0[indx1][ii*3+jj]=0;
					}
				break;

			default:
				break;
			}
		}
	return 0;
}

int Gelastic::Set_SFTS_habit()
{
//	cm=habit plane normal
//	cb=average transformation direction
//	transformationstrain_ij=0.1308*0.5*(b_i*m_j+b_j*m_i)
//	parameter(nv=24)
/*
	float a1,a2,a3,a4,a5,a6;
	float cm[4][25]={0},cb[4][25]={0};
	int i;
	int ii,jj;
	long int indx1;

	if(nv!=24)
	{
		cout<<"Invalid number of variants for habit-v!"<<endl;
		exit(1);
	}

	a1=0.8888;
	a2=0.4045;
	a3=0.2153;
	a4=0.4343;
	a5=0.4878;
	a6=0.7576;
	cm[1][1]=-a1;
	cm[2][1]=-a2;
	cm[3][1]=a3;
	cb[1][1]=a4;
	cb[2][1]=-a5;
	cb[3][1]=a6;
      
//	m[1]
	for(i=1;i<=8;i++)
         cm[1][i]=-a1;
	cm[1][9]=a2;
	cm[1][10]=-a2;
	cm[1][11]=a3;
	cm[1][12]=a3;
	
	cm[1][13]=-a3;
	cm[1][14]=-a3;
	cm[1][15]=a2;
	cm[1][16]=-a2;
	cm[1][17]=a3;
	cm[1][18]=a3;
	cm[1][19]=a2;
	cm[1][20]=-a2;
	cm[1][21]=a2;
	cm[1][22]=-a2;
	cm[1][23]=-a3;
	cm[1][24]=-a3;

//	m[2]
	cm[2][1]=-a2;
	cm[2][2]=a2;
	cm[2][3]=a3;
	cm[2][4]=a3;
	cm[2][5]=-a3;
	cm[2][6]=-a3;
	cm[2][7]=a2;
	cm[2][8]=-a2;
	cm[2][9]=-a1;
	cm[2][10]=-a1;
	cm[2][11]=-a1;
	cm[2][12]=-a1;

	cm[2][13]=-a1;
	cm[2][14]=-a1;
	cm[2][15]=-a1;
	cm[2][16]=-a1;
	cm[2][17]=a2;
	cm[2][18]=-a2;
	cm[2][19]=a3;
	cm[2][20]=a3;
	cm[2][21]=-a3;
	cm[2][22]=-a3;
	cm[2][23]=a2;
	cm[2][24]=-a2;


//	m[3]
	cm[3][1]=a3;
	cm[3][2]=a3;
	cm[3][3]=-a2;
	cm[3][4]=a2;
	cm[3][5]=a2;
	cm[3][6]=-a2;
	cm[3][7]=-a3;
	cm[3][8]=-a3;
	cm[3][9]=a3;
	cm[3][10]=a3;
	cm[3][11]=-a2;
	cm[3][12]=a2;


	cm[3][13]=a2;
	cm[3][14]=-a2;
	cm[3][15]=-a3;
	cm[3][16]=-a3;
	for(i=17;i<=24;i++)
		cm[3][i]=-a1;
 
//	b[1]
	for(i=1;i<=8;i++)
		cb[1][i]=a4;
	cb[1][9]=a5;
	cb[1][10]=-a5;
	cb[1][11]=a6;
	cb[1][12]=a6;

	cb[1][13]=-a6;
	cb[1][14]=-a6;
	cb[1][15]=a5;
	cb[1][16]=-a5;
	cb[1][17]=a6;
	cb[1][18]=a6;
	cb[1][19]=a5;
	cb[1][20]=-a5;
	cb[1][21]=a5;
	cb[1][22]=-a5;
	cb[1][23]=-a6;
	cb[1][24]=-a6;
      

//	b[2]
	cb[2][1]=-a5;
	cb[2][2]=a5;
	cb[2][3]=a6;
	cb[2][4]=a6;
	cb[2][5]=-a6;
	cb[2][6]=-a6;
	cb[2][7]=a5;
	cb[2][8]=-a5;
	for(i=9;i<=16;i++)
		cb[2][i]=a4;
		
	cb[2][17]=a5;
	cb[2][18]=-a5;
	cb[2][19]=a6;
	cb[2][20]=a6;
	cb[2][21]=-a6;
	cb[2][22]=-a6;
	cb[2][23]=a5;
	cb[2][24]=-a5;
	
//	b[3]
	cb[3][1]=a6;
	cb[3][2]=a6;
	cb[3][3]=-a5;
	cb[3][4]=a5;
	cb[3][5]=a5;
	cb[3][6]=-a5;
	cb[3][7]=-a6;
	cb[3][8]=-a6;
	cb[3][9]=a6;
	cb[3][10]=a6;
	cb[3][11]=-a5;
	cb[3][12]=a5;

	cb[3][13]=a5;
	cb[3][14]=-a5;
	cb[3][15]=-a6;
	cb[3][16]=-a6;
	for(i=17;i<=24;i++)
		cb[3][i]=a4;

	for(i=0;i<24;i++)
	{
		indx1=index1(0,i);
		for(ii=0;ii<3;ii++)
			for(jj=0;jj<3;jj++)
			{
				(e0[indx1])[ii*3+jj]=0.1308*0.5*(cb[ii+1][i+1]*cm[jj+1][i+1]+cb[jj+1][i+1]*cm[ii+1][i+1]);
			}
	}
*/
	return 0;
}

void Gelastic::LU_dcmp(real **a, int n, int *indx, real *d)
{

	/* Adopted from "Numerical Recipeis in C" by W Press et. al.

	   Given a matrix a[1..n][1..n], this function REPLACE it by
	   the LU decomposition of a rowwise permutation of itself.
INPUT:
	a -- original matrix
	n -- dimension
OUTPUT:
	a -- matrix containing L and U
	indx -- row permuation effected by the partial pivoting
	d -- +1/-1, depending on the number of row interchanges was even or odd, respectively. */

	int i,imax,j,k;
	real big,dum,sum,temp;
	real *vv; //vv stores the implicit scaling of each row.

	vv=(real*)malloc(n*sizeof(real));
	*d=1.0;		//No row interchanges yet.
	for (i=1;i<=n;i++) { //Loop over rows to get the implicit scaling information
		big=0.0; 
		for (j=1;j<=n;j++)
		if ((temp=fabs(a[i-1][j-1])) > big) big=temp;
		if (fabs(big)<1E-10){
			printf("Singular matrix in function LU_dcmp()!!\n");
			for(j=1;j<=n;j++){
					printf("LU_dcmp: a[%d][%d]=%lf\n",i,j,a[i-1][j-1]);
					fflush(stdout);
			}
			exit(35);
		}
		//No nonzero largest element.
		vv[i-1]=1.0/big; //Save the scaling.
	}
	for (j=1;j<=n;j++) {	//This is the loop over columns of Crout's method.
		for (i=1;i<j;i++) {		//This is equation (2.3.12) except for i = j.
			sum=a[i-1][j-1];
			for (k=1;k<i;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
		}
		big=0.0;	//Initialize for the search for largest pivot element.
		for (i=j;i<=n;i++) {	//This is i = j of equation (2.3.12) and i = j+1...N
			sum=a[i-1][j-1];		//of equation (2.3.13).
			for (k=1;k<j;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
			if ( (dum=vv[i-1]*fabs(sum)) >= big) {
				// Is the figure of merit for the pivot better than the best so far?
				big=dum;
				imax=i;
			}
		}
		//Do we need to interchange rows?
		if (j != imax) {	//Yes, do so...
			for (k=1;k<=n;k++) {	
				dum=a[imax-1][k-1];
				a[imax-1][k-1]=a[j-1][k-1];
				a[j-1][k-1]=dum;
			}
			*d = -(*d);		//...and change the parity of d.
			vv[imax-1]=vv[j-1];	//Also interchange the scale factor.
		}
		indx[j-1]=imax;
		if (fabs(a[j-1][j-1]) < 1E-5) a[j-1][j-1]=TINY;
//		 If the pivot element is zero the matrix is singular (at least to the precision of the
//		algorithm). For some applications on singular matrices, it is desirable to substitute
//		TINY for zero.
		if (j != n) {	//Now, finally, divide by the pivot element.
			dum=1.0/(a[j-1][j-1]);
			for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
		}
	}	//Go back for the next column in the reduction.
	
	free(vv);

	return;
}
/*end LU_dcmp()*/

void Gelastic::LU_bksb(real **a, int n, int *indx, real b[])
{
	/* Adopted from "Numerical Recipeis in C" by W Press et. al.

	   Solve the set of n linear equations A*X = B.
INPUT:
	a -- matrix A in the LU decomposition form obtained through LU_dcmp()
	n -- dimension
	indx -- output of LU_dcmp()
	b -- vextor B
OUTPUT:
	b -- solution X.							*/

	int i,ii=0,ip,j;
	real sum;

//	 When ii is set to a positive value, it will become the
//		index of the first nonvanishing element of b. We now
//		do the forward substitution, equation (2.3.6). The
//		only new wrinkle is to unscramble the permutation
//		as we go.
	for (i=1;i<=n;i++) {
		ip=indx[i-1];
		sum=b[ip-1];
		b[ip-1]=b[i-1];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
		else if (fabs(sum)>1E-5) ii=i;		//A nonzero element was encountered, so from now on we
								//will have to b[i]=sum; do the sums in the loop above.
		b[i-1] = sum;
	}
	for (i=n;i>=1;i--) {	//Now we do the backsubstitution, equation (2.3.7).
		sum=b[i-1];
		for (j=i+1;j<=n;j++) sum -= a[i-1][j-1]*b[j-1];
		b[i-1]=sum/a[i-1][i-1];	//Store a component of the solution vector X.
	}	// All done!
}
/*end LU_bksb()*/

void Gelastic::LU_inv_66(voigt66 c)
{
//	 Inverse the matrix using LU decomposition.
//	   The original matrix will be replaced with
//	   its inverse.
//		*A special case for 6x6 matrix with voigt66 type,
//		modified from LU_inverse(real**, int) 

	int indx[6];
	int n;
	real d;
	int i, j;
	voigt col;
	voigt66 y;
	real *a[6] = {c[0],c[1],c[2],c[3],c[4],c[5]};

	n = 6;
	LU_dcmp(a,n,indx,&d);
	for(j=0;j<n;j++){
		for(i=0;i<n;i++) col[i] = 0.0;
		col[j] = 1.0;
		LU_bksb(a,n,indx,col);
		for(i=0;i<n;i++) y[i][j] = col[i];
	}

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			c[i][j] = y[i][j];

	return;
}

/*end LU_inv_66()*/
/*
static inline void Gelastic::chg_basis_Kelvin_4(ten4th T4, voigt66 C2)
    {
       C4_loop{
            int p = (mi+1)*(mi==mj)+(1-(mi==mj))*(7-mi-mj);
            int q = (mk+1)*(mk==ml)+(1-(mk==ml))*(7-mk-ml);
            real t1 = ((real)(mi==mj)+(1.-(mi==mj))/RSQ2);
            real t2 = ((real)(mk==ml)+(1.-(mk==ml))/RSQ2);
            C2[p-1][q-1] = t1*t2*T4[mi][mj][mk][ml];
        }
    }

static inline void Gelastic::chg_basis_Kelvin_3( voigt66 C2, ten4th T4)
{
        C4_loop{
            int p = (mi+1)*(mi==mj)+(1-(mi==mj))*(7-mi-mj);
            int q = (mk+1)*(mk==ml)+(1-(mk==ml))*(7-mk-ml);
            real t1 = ((real)(mi==mj)+(1.-(mi==mj))/RSQ2);
            real t2 = ((real)(mk==ml)+(1.-(mk==ml))/RSQ2);
            T4[mi][mj][mk][ml] = C2[p-1][q-1]/t1/t2;
        }
}
*/
