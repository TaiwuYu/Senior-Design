//////////////////////////
// Ggradient.cpp //////
// Y.Gao 05-18-2009	//////
//////////////////////////

#include "Ggradient.h"

Ggradient::Ggradient()
{
	ng=nv=nx=ny=nz=0;
	kappa=NULL;
	initflag=false;
	option=ISO_P7;
}

Ggradient::~Ggradient()
{
	if(initflag==true)
		delete[] kappa;
	kappa=NULL;
}

int Ggradient::Set(int n[],float *ka,gradient_option op)
{
	Set_n(n);
	Set_kappa(ka);
	option=op;
	initflag=true;

	return 0;
}

int Ggradient::Set_n(int n[])
{
	if(n[0]<1 || n[1]<1 || n[2]<1 || n[3]<1 || n[4]<1)
	{
		cout<<"Invalid dimensions in Ggradient!"<<endl;
		exit(1);
	}
	ng=n[0];
	nv=n[1];
	nx=n[2];
	ny=n[3];
	nz=n[4];

	return 0;
}

int Ggradient::Set_kappa(float *ka)
{
	int i;

	kappa=new float[nv*nv];
	assert(kappa!=NULL);

	for(i=0;i<nv*nv;i++)
		kappa[i]=ka[i];

	return 0;
}

float Ggradient::Energy(float* eta[],int *gs)
{
	float result=0;
	switch(option)
	{
	case ISO_P7:
		result=E_ISO_P7(eta,gs);
		break;
	case ISO_P27:
		cout<<"No ISO_P27 option in Ggradient!"<<endl;
		exit(1);
//		result=E_ISO_P27(eta,gs);
//		break;
	case ANI_P7:
		cout<<"No ANI_P7 option in Ggradient!"<<endl;
		exit(1);
//		result=E_ANI_P7(eta,gs);
//		break;
	case ANI_P27:
		cout<<"No ANI_P27 option in Ggradient!"<<endl;
		exit(1);
//		result=E_ANI_P27(eta,gs);
//		break;
	default:
		cout<<"Invalid Ggradient option!"<<endl;
		exit(1);
	}
	return result;
}

int Ggradient::Potential(float* eta[],float* d_eta[],int *gs)
{
	switch(option)
	{
	case ISO_P7:
		Miu_ISO_P7(eta,d_eta,gs);
		break;
	case ISO_P27:
		Miu_ISO_P27(eta,d_eta,gs);
		break;
	case ANI_P7:
		Miu_ANI_P7(eta,d_eta,gs);
		break;
	case ANI_P27:
		Miu_ANI_P27(eta,d_eta,gs);
		break;
	default:
		cout<<"Invalid Ggradient option!"<<endl;
		exit(1);
	}
	return 0;
}

float Ggradient::E_ISO_P7(float* eta[],int *gs)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	float xp,xm,yp,ym,zp,zm;
	float kap=100*kappa[0];

	float result=0;
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
						
						indx2=index2((i+1)%nx,j,k);
						xp=(eta[indx1])[indx2];
						indx2=index2((i-1+nx)%nx,j,k);
						xm=(eta[indx1])[indx2];

						indx2=index2(i,(j+1)%ny,k);
						yp=(eta[indx1])[indx2];
						indx2=index2(i,(j-1+ny)%ny,k);
						ym=(eta[indx1])[indx2];

						indx2=index2(i,j,(k+1)%nz);
						zp=(eta[indx1])[indx2];
						indx2=index2(i,j,(k-1+nz)%nz);
						zm=(eta[indx1])[indx2];
						
						switch(v1)
						{
						case 0:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 1:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 2:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 3:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 6:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 7:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 4:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						case 5:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
                                                case 8:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
                                                case 9:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
                                                case 10:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
                                                case 11:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
                                                case 12:
							result+=kappa[0]*((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
							break;
						default:
							cout<<"Invalid variant ID in Ggradient E_ISO_P7!"<<endl;
							break;
						}
					}
				}
			}
	result/=8.0;
	return result;
}

float Ggradient::E_ISO_P27(float* eta[],int *gs)
{
	float result=0;
	return result;
}

float Ggradient::E_ANI_P7(float* eta[],int *gs)
{
	float result=0;
	return result;
}

float Ggradient::E_ANI_P27(float* eta[],int *gs)
{
	float result=0;
	return result;
}
int Ggradient::Miu_ISO_P7(float* eta[],float* d_eta[],int *gs)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	float xp,xm,yp,ym,zp,zm,temp;

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
						indx2=index2(i,j,k);
						temp=eta[indx1][indx2];


						indx2=index2((i+1)%nx,j,k);
						xp=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,j,k);
						xm=(eta[indx1])[indx2];

						indx2=index2(i,(j+1)%ny,k);
						yp=(eta[indx1])[indx2];
						indx2=index2(i,(j-1+ny)%ny,k);
						ym=(eta[indx1])[indx2];

						indx2=index2(i,j,(k+1)%nz);
						zp=(eta[indx1])[indx2];
						indx2=index2(i,j,(k-1+nz)%nz);
						zm=(eta[indx1])[indx2];
						
						indx2=index2(i,j,k);
                             //                   if(v1>11)
			//			continue;
				//		if(v1>11)
				//			d_eta[v1][indx2]-=0.001*kappa[0]*(xp+xm+yp+ym+zp+zm-6*temp);
/*				if(v1==6)
							d_eta[v1][indx2]-=10*kappa[0]*(xp+xm+yp+ym+zp+zm-6*temp);
					else if(v1==7)
							d_eta[v1][indx2]-=5*kappa[0]*(xp+xm+yp+ym+zp+zm-6*temp);
                                       else{
				       */
							d_eta[v1][indx2]-=kappa[0]*(xp+xm+yp+ym+zp+zm-6*temp);
				//	}
					}
				}
			}

	return 0;
}

int Ggradient::Miu_ISO_P27(float* eta[],float* d_eta[],int *gs)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	float temp,
		xp,xm,yp,ym,zp,zm,tp1,
		ppo,pmo,mpo,mmo,pop,pom,mop,mom,opp,opm,omp,omm,tp2,
		ppp,ppm,pmp,pmm,mpp,mpm,mmp,mmm,tp3;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
					for(v1=0;v1<nv;v1++)	// Only in single crystal
					{
						indx1=index1(g1,v1);
						indx2=index2(i,j,k);
						temp=(eta[indx1])[indx2];

//	Nearest Neighbor
						indx2=index2((i+1)%nx,j,k);
						xp=(eta[indx1])[indx2];
						indx2=index2((i-1+nx)%nx,j,k);
						xm=(eta[indx1])[indx2];

						indx2=index2(i,(j+1)%ny,k);
						yp=(eta[indx1])[indx2];
						indx2=index2(i,(j-1+ny)%ny,k);
						ym=(eta[indx1])[indx2];

						indx2=index2(i,j,(k+1)%nz);
						zp=(eta[indx1])[indx2];
						indx2=index2(i,j,(k-1+nz)%nz);
						zm=(eta[indx1])[indx2];

						tp1=xp+xm+yp+ym+zp+zm;
//	Second Nearest Neighbor
						indx2=index2((i+1)%nx,(j+1)%ny,k);
						ppo=eta[indx1][indx2];
						indx2=index2((i+1)%nx,(j-1+ny)%ny,k);
						pmo=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,(j+1)%ny,k);
						mpo=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,(j-1+ny)%ny,k);
						mmo=eta[indx1][indx2];

						indx2=index2((i+1)%nx,j,(k+1)%nz);
						pop=eta[indx1][indx2];
						indx2=index2((i+1)%nx,j,(k-1+nz)%nz);
						pom=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,j,(k+1)%nz);
						mop=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,j,(k-1+nz)%nz);
						mom=eta[indx1][indx2];

						indx2=index2(i,(j+1)%ny,(k+1)%nz);
						opp=eta[indx1][indx2];
						indx2=index2(i,(j+1)%ny,(k-1+nz)%nz);
						opm=eta[indx1][indx2];
						indx2=index2(i,(j-1+ny)%ny,(k+1)%nz);
						omp=eta[indx1][indx2];
						indx2=index2(i,(j-1+ny)%ny,(k-1+nz)%nz);
						omm=eta[indx1][indx2];
						
						tp2=ppo+pmo+mpo+mmo+pop+pom+mop+mom+opp+opm+omp+omm;
//	Third Nearest Neighbor
						indx2=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
						ppp=eta[indx1][indx2];
						indx2=index2((i+1)%nx,(j+1)%ny,(k-1+nz)%nz);
						ppm=eta[indx1][indx2];

						indx2=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
						pmp=eta[indx1][indx2];
						indx2=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
						pmm=eta[indx1][indx2];

						indx2=index2((i-1+nx)%nx,(j+1)%ny,(k+1)%nz);
						mpp=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
						mpm=eta[indx1][indx2];

						indx2=index2((i-1+nx)%nx,(j-1+ny)%ny,(k+1)%nz);
						mmp=eta[indx1][indx2];
						indx2=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
						mmm=eta[indx1][indx2];

						tp3=ppp+ppm+pmp+pmm+mpp+mpm+mmp+mmm;
						
						indx2=index2(i,j,k);
						(d_eta[v1])[indx2]-=kappa[0]*(0.5*tp1+0.125*tp2-4.5*temp);
					}
				}
			}

	return 0;
}

int Ggradient::Miu_ANI_P7(float* eta[],float* d_eta[],int *gs)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	float xp,xm,yp,ym,zp,zm,temp,tp0,tp1,tp2;
	float coef=1;

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
						indx2=index2(i,j,k);
						temp=(eta[indx1])[indx2];


						indx2=index2((i+1)%nx,j,k);
						if(gs[indx2]!=g1)
							xp=temp;
						else
							xp=(eta[indx1])[indx2];
						indx2=index2((i-1+nx)%nx,j,k);
						if(gs[indx2]!=g1)
							xm=temp;
						else
							xm=(eta[indx1])[indx2];

						indx2=index2(i,(j+1)%ny,k);
						if(gs[indx2]!=g1)
							yp=temp;
						else
							yp=(eta[indx1])[indx2];
						indx2=index2(i,(j-1+ny)%ny,k);
						if(gs[indx2]!=g1)
							ym=temp;
						else
							ym=(eta[indx1])[indx2];

						indx2=index2(i,j,(k+1)%nz);
						if(gs[indx2]!=g1)
							zp=temp;
						else
							zp=(eta[indx1])[indx2];
						indx2=index2(i,j,(k-1+nz)%nz);
						if(gs[indx2]!=g1)
							zm=temp;
						else
							zm=(eta[indx1])[indx2];

						
						indx2=index2(i,j,k);
// modified for no penalty on y-z plane
						(d_eta[v1])[indx2]-=kappa[0]*(yp+ym+zp+zm-4*temp);
//	Among 0 1 2 Variants
						switch(v1)
						{
						case 0:
							tp0=kappa[0]*(xp+xm-2*temp);
							break;
						case 1: 
							tp1=kappa[0]*(xp+xm-2*temp);
							break;
						case 2:
							tp2=kappa[0]*(xp+xm-2*temp);
							break;
						default:
							break;
						}
					}
					d_eta[0][indx2]-=coef*(tp1+tp2);
					d_eta[1][indx2]-=coef*(tp0+tp2);
					d_eta[2][indx2]-=coef*(tp0+tp1);
				}
			}

	return 0;
}

int Ggradient::Miu_ANI_P27(float* eta[],float* d_eta[],int *gs)
{
	int g1,v1,g2,v2,i,j,k;
	int ii,jj;
	long int indx1_1,indx1_2,indx2,indx22;
	float	e1,												// itself
			poo,moo,opo,omo,oop,oom,						// 1st nearest
			ppo,pmo,mpo,mmo,pop,pom,mop,mom,opp,opm,omp,omm,// 2nd nearest
			ppp,ppm,pmp,pmm,mpp,mpm,mmp,mmm;				// 3rd nearest
	float temp;
	float coef=0.995;

	float kap[16][16][3][3]={0}; // maximum ng*nv<=16

//	For lower energy in {111}
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1_1=index1(g1,v1);
			for(g2=0;g2<ng;g2++)
				for(v2=0;v2<nv;v2++)
				{
					indx1_2=index1(g2,v2);
					for(ii=0;ii<3;ii++)
						for(jj=0;jj<3;jj++)
						{
							if(indx1_1==indx1_2 && ii==jj)
								kap[indx1_1][indx1_2][ii][jj]=kappa[0];
/*							if(indx1_1!=indx1_2)
							{
								if(indx1_1<3 && indx1_2<3 && indx1_1>=0 && indx1_2>=0) //(111)
									kap[indx1_1][indx1_2][ii][jj]=kappa[0]/3.0*coef;
								if(indx1_1<6 && indx1_2<6 && indx1_1>=3 && indx1_2>=3) //(-111)
								{
									if( (ii==0 && jj!=0) || (ii!=0 && jj==0) )
										kap[indx1_1][indx1_2][ii][jj]=-kappa[0]/3.0*coef;
									else
										kap[indx1_1][indx1_2][ii][jj]=kappa[0]/3.0*coef;
								}
								if(indx1_1<9 && indx1_2<9 && indx1_1>=6 && indx1_2>=6) //(1-11)
								{
									if( (ii==1 && jj!=1) || (ii!=1 && jj==1) )
										kap[indx1_1][indx1_2][ii][jj]=-kappa[0]/3.0*coef;
									else
										kap[indx1_1][indx1_2][ii][jj]=kappa[0]/3.0*coef;
								}
								if(indx1_1<12 && indx1_2<12 && indx1_1>=9 && indx1_2>=9) //(11-1)
								{
									if( (ii==2 && jj!=2) || (ii!=2 && jj==2) )
										kap[indx1_1][indx1_2][ii][jj]=-kappa[0]/3.0*coef;
									else
										kap[indx1_1][indx1_2][ii][jj]=kappa[0]/3.0*coef;
								}

//								if(ii<2 && jj<2) // used for (110) test
//									kap[indx1_1][indx1_2][ii][jj]=kappa[0]/2.0*coef;

							}
*/						}
					if(indx1_1!=indx1_2) // used for (100) test
						kap[indx1_1][indx1_2][0][0]=kappa[0]*coef;	
				}
		}

//	Gradient potential Calculation
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
					for(v1=0;v1<nv;v1++) // eta[s]: g1,v1
					{
						indx1_1=index1(g1,v1);
						temp=0;
						for(g2=0;g2<ng;g2++)
						{
							if(g2!=gs[indx2])
								continue;
							for(v2=0;v2<nv;v2++) // eta[p]: g2,v2
							{
								indx1_2=index1(g2,v2);
								e1=eta[indx1_2][indx2];
								for(ii=0;ii<3;ii++)
									for(jj=0;jj<3;jj++)
									{
									//	2nd-order partial derivative
										if(ii==0)
										{
											switch(jj)
											{
											case 0:	//xx
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2((i+1)%nx,j,k);
												poo=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,j,k);
												moo=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]*(poo+moo-2*e1);
												break;

											case 1: //xy
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2((i+1)%nx,(j+1)%ny,k);
												ppo=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,k);
												mmo=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,k);
												pmo=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,k);
												mpo=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/8.0*(ppo+mmo-pmo-mpo);


												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k+1)%nz);
												mmp=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
												pmp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k+1)%nz);
												mpp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppp+mmp-pmp-mpp);

												
												indx22=index2((i+1)%nx,(j+1)%ny,(k-1+nz)%nz);
												ppm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												pmm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mpm=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppm+mmm-pmm-mpm);

												break;

											case 2: //xz
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2((i+1)%nx,j,(k+1)%nz);
												pop=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,j,(k-1+nz)%nz);
												mom=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,j,(k+1)%nz);
												pom=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,j,(k-1+nz)%nz);
												mop=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/8.0*(pop+mom-pom-mop);


												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
												mpm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
												mpp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppp+mpm-ppm-mpp);


												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
												pmp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
												pmm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(pmp+mmm-pmm-mmp);
												
												break;

											default:
												break;
											}
										}
										if(ii==1)
										{
											switch(jj)
											{
											case 0: //yx
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2((i+1)%nx,(j+1)%ny,k);
												ppo=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,k);
												mmo=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,k);
												pmo=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,k);
												mpo=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/8.0*(ppo+mmo-pmo-mpo);


												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k+1)%nz);
												mmp=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
												pmp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k+1)%nz);
												mpp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppp+mmp-pmp-mpp);

												
												indx22=index2((i+1)%nx,(j+1)%ny,(k-1+nz)%nz);
												ppm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												pmm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mpm=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppm+mmm-pmm-mpm);

												break;

											case 1: //yy
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2(i,(j+1)%ny,k);
												opo=eta[indx1_2][indx22];

												indx22=index2(i,(j-1+ny)%ny,k);
												omo=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]*(opo+omo-2*e1);
												break;

											case 2: //yz
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2(i,(j+1)%ny,(k+1)%nz);
												opp=eta[indx1_2][indx22];

												indx22=index2(i,(j-1+ny)%ny,(k-1+nz)%nz);
												omm=eta[indx1_2][indx22];

												indx22=index2(i,(j+1)%ny,(k-1+nz)%nz);
												opm=eta[indx1_2][indx22];

												indx22=index2(i,(j-1+ny)%ny,(k-1+nz)%nz);
												omp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/8.0*(opp+omm-opm-omp);


												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppp=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												pmm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j+1)%ny,(k-1+nz)%nz);
												ppm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												pmp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppp+pmm-ppm-pmp);


												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k+1)%nz);
												mpp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
												mpm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(mpp+mmm-mpm-mmp);

												break;

											default:
												break;
											}
										}
										if(ii==2)
										{
											switch(jj)
											{
											case 0: //zx
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2((i+1)%nx,j,(k+1)%nz);
												pop=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,j,(k-1+nz)%nz);
												mom=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,j,(k+1)%nz);
												pom=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,j,(k-1+nz)%nz);
												mop=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/8.0*(pop+mom-pom-mop);


												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
												mpm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
												mpp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppp+mpm-ppm-mpp);


												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
												pmp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k+1)%nz);
												pmm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(pmp+mmm-pmm-mmp);

												break;

											case 1: //zy
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2(i,(j+1)%ny,(k+1)%nz);
												opp=eta[indx1_2][indx22];

												indx22=index2(i,(j-1+ny)%ny,(k-1+nz)%nz);
												omm=eta[indx1_2][indx22];

												indx22=index2(i,(j+1)%ny,(k-1+nz)%nz);
												opm=eta[indx1_2][indx22];

												indx22=index2(i,(j-1+ny)%ny,(k-1+nz)%nz);
												omp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/8.0*(opp+omm-opm-omp);


												indx22=index2((i+1)%nx,(j+1)%ny,(k+1)%nz);
												ppp=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												pmm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j+1)%ny,(k-1+nz)%nz);
												ppm=eta[indx1_2][indx22];

												indx22=index2((i+1)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												pmp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(ppp+pmm-ppm-pmp);


												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k+1)%nz);
												mpp=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j+1)%ny,(k-1+nz)%nz);
												mpm=eta[indx1_2][indx22];

												indx22=index2((i-1+nx)%nx,(j-1+ny)%ny,(k-1+nz)%nz);
												mmp=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]/16.0*(mpp+mmm-mpm-mmp);

												break;

											case 2: //zz
												if(kap[indx1_1][indx1_2][ii][jj]<0.00001 && kap[indx1_1][indx1_2][ii][jj]>-0.00001)
													break;
												indx22=index2(i,j,(k+1)%nz);
												oop=eta[indx1_2][indx22];

												indx22=index2(i,j,(k-1+nz)%nz);
												oom=eta[indx1_2][indx22];

												temp-=kap[indx1_1][indx1_2][ii][jj]*(oop+oom-2*e1);
												break;

											default:
												break;
											}
										}
									}
							}
						}

						(d_eta[v1])[indx2]+=temp;
					}
				}
			}

	return 0;
} 

