////////////////////////
// GInitConfig.cpp	////
// Y.Gao 05-20-2009	////
////////////////////////

#include "GInitConfig.h"
#include "shortcut_wrappers.h"
#include "Gmt.h"
#include <iomanip>
#include <string.h>
#include <stdlib.h>

GInitConfig::GInitConfig()
{
	ng=nv=nx=ny=nz=0;
	initflag=false;
	option=HOMO;
}


GInitConfig::~GInitConfig()
{

}

int GInitConfig::Set_n(int n[])
{
	if(n[0]<1 || n[1]<1 || n[2]<1 || n[3]<1 || n[4]<1)
	{
		cout<<"Invalid dimensions in GInitConfig!"<<endl;
		1;
	}
	ng=n[0];
	nv=n[1];
	nx=n[2];
	ny=n[3];
	nz=n[4];

	return 0;
}

int GInitConfig::Set(int n[],init_option op,
					 float *eta[],float a[],char eta_file[],
					 int *gs,char gs_file[],float *conc1,float *conc2)
{
	Set_n(n);
	option=op;
long int indx1,indx2;
int g1,v1,i,j,k;
	Init_gs(gs,gs_file);
	switch(option)
	{
	case HOMO:
		Init_HOMO(eta);
		break;
	case SINGLE:
		Init_SINGLE(eta,a,gs);
		break;
	case TWIN_13:
		Init_TWIN_13(eta,a,gs);
		break;
	case USER:
		Init_USER(eta,a,gs,conc1,conc2);
		break;
	case LOAD:
		Init_LOAD(eta,eta_file);
		break;
	case GRAIN_3:
		Init_GRAIN_3(eta);
		break;
	case GRAIN_4:
		Init_GRAIN_4(eta);
		break;
	case DOUBLE:
		Init_DOUBLE(eta,a,gs);
		break;
	case TWIN_12:
		Init_TWIN_12(eta,a,gs);
		break;
	case TWIN_159:
		Init_TWIN_159(eta,a,gs);
		break;
	case LAM_123:
		Init_LAM_123(eta,eta_file,a,gs);
		break;
	case LOAD_3264C:
		Init_LOAD_3264C(eta,eta_file);
		break;
	case TWO:
		Init_TWO(eta,a,gs);
		break;
	case HPHASE:
		Init_Hphase(eta,a,gs);
		break;
        case USER_SET:
             //   Init_field(eta,input_eta);
                Init_field(eta,conc1,conc2);
                break;
        case RANDOM_PARTICLE:
	        Init_particle(eta,a,gs);
		break;

	default:
		cout<<"Invalid GInitConfig option!"<<endl;
		1;
	}
	
						for(g1=0;g1<ng;g1++)
						for(v1=0;v1<nv;v1++)
						{
               if(v1<10)
               {
	                          char fn[]="eta_init_0.vtk";
							indx1=index1(g1,v1);
			fn[9]=(char)((int)'0'+g1*nv+v1);
			ofstream fout(fn,ios::out);
			Output_VTK_header(&fout,nx,ny,nz);
			fout<<fixed;
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						fout<<setprecision(6)<<(eta[indx1])[indx2]<<endl;
						//fout<<(eta[indx1])[indx2]<<endl;
					}
                }
             else{
	                          char fn[]="eta_init_10.vtk";
							indx1=index1(g1,v1);
			fn[10]=(char)((int)'0'+g1*nv+v1-10);
			ofstream fout(fn,ios::out);
			Output_VTK_header(&fout,nx,ny,nz);
			fout<<fixed;
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						fout<<setprecision(6)<<(eta[indx1])[indx2]<<endl;
					}
                          
                      }
		}
		
	initflag=true;
	return 0;
}

int GInitConfig::Init_HOMO(float *eta[])	// Ni30Pt20Ti50
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					for(k=0;k<nz;k++)
					{
						indx2=index2(i,j,k);
						switch(v1)
						{
						case 0:
							eta[indx1][indx2]=0;
							break;
						case 1:
							eta[indx1][indx2]=0;
							break;
						case 2:
							eta[indx1][indx2]=0;
							break;
						case 3:
							eta[indx1][indx2]=0;
							break;
						case 4:
							eta[indx1][indx2]=0;
							break;
						case 5:
							eta[indx1][indx2]=0;
							break;
						case 6://X_Ni
					//		eta[indx1][indx2]=0.57;
							eta[indx1][indx2]=0.0;
							break;
						case 7://X_Hf
						//	eta[indx1][indx2]=0.27;
							eta[indx1][indx2]=0.0;
							break;
						case 8:
						case 9:
						case 10:
						case 11:
						case 12:
						case 13:
						case 14:
						case 15:
						case 16:
						case 17:
							eta[indx1][indx2]=0;
							break;
			//			case 4:
			//				eta[indx1][indx2]=0;
						//	eta[indx1][indx2]=0;
			//				break;
			//			case 5:
			//				eta[indx1][indx2]=0;
						//	eta[indx1][indx2]=0;
			//				break;
						default:
							cout<<"Invalid variant ID!_Config"<<endl;
							//exit(1);
						}
					}
		}
	return 0;
}

int GInitConfig::Init_SINGLE(float *eta[],float a[],int *gs) // 1D test for martensite
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	float init_vol=4.0/3.0*3.1416*a[0]*a[0]*a[0];

	if(nz==1)	// 2D
	{
		init_vol=3.1416*a[0]*a[0];
		if(ny==1)	// 1D
			init_vol=2*a[0]+1;
	}

	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1==gs[indx2])
					{
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
                                                   
							switch(v1)
							{
					    case 0:
						//		if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
						//			(eta[indx1])[indx2]=1;
								break;
					    case 1:
								break;
					    case 2:
								break;
					    case 3:
								break;
					    case 4:
								if(abs(i-nx/2)<=a[0]/2||i<=a[0]/2||(nx-i)<=a[0]/2||k<=a[0]/2||(nz-k)<=a[0]/2||j<=a[0]/2||(ny-j)<=a[0]/2)
								//if(abs(i-nx/2)<=a[0]/2||i<=a[0]/2||(nx-i)<=a[0]/2||k<=a[0]/2||(nz-k)<=a[0]/2)
									(eta[indx1])[indx2]=1;
								break;
					    case 5:
								break;
					    case 6:
								break;
					    case 7:
								break;
					    case 8:
								break;
					    case 9:
								break;
					    case 10:
								break;
					    case 11:
								break;
					    case 12:
								if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
								//if((i-nx/2)<0 && (i-nx/4)>0)
						//		&& fabs((float)(i-nx/2)*(0.707)+(float)(j-ny/2)*(-0.707))<=a[1])
						//		if(fabs(i-nx/2)<a[0]/2 && fabs(j-ny/2)<a[0]/2 && fabs(k-nz/2)<a[0]/2)
									(eta[indx1])[indx2]=0;
						//		else if((i-3*nx/4)>0)
						//			(eta[1])[indx2]=1.0;
                                                //                   else
						//			(eta[indx1])[indx2]=0;
								break;
	//						case 4:
							/*	if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
								//&& fabs((float)(i-nx/2)*(0.707)+(float)(j-ny/2)*(-0.707))<=a[1])
//								if(fabs(i-nx/2)<a[0]/2 && fabs(j-ny/2)<a[0]/2 && fabs(k-nz/2)<a[0]/2)
									eta[indx1][indx2]=0.375;
								else
									eta[indx1][indx2]=(nx*ny*nz*0.3-init_vol*0.375)/(nx*ny*nz-init_vol);
							*/
          //                                                        	break;
	//						case 5:
							/*	if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
								// && fabs((float)(i-nx/2)*(0.707)+(float)(j-ny/2)*(-0.707))<=a[1])
//								if(fabs(i-nx/2)<a[0]/2 && fabs(j-ny/2)<a[0]/2 && fabs(k-nz/2)<a[0]/2)
									eta[indx1][indx2]=0.167;
								else
									eta[indx1][indx2]=(nx*ny*nz*0.2-init_vol*0.167)/(nx*ny*nz-init_vol);
								*/
          //                                                               break;

							default:
								(eta[indx1])[indx2]=0;
							}//end of variant
					}
    }
					else
					{
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
							(eta[indx1])[indx2]=0;
						}
					}
				}
			}//end of space
		
	return 0;
}

int GInitConfig::Init_TWO(float *eta[],float a[],int *gs)
{
	int g1,v1,i,j,k,ix,iy,iz,jx,jy,jz;
	long int indx1,indx2;
	int iv,ih;
	float tp1,tp2,tp3,mp1,mp2,mp3;

	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
		iv=0;
	//iv=0;
//		ih=rand()%2;
		ih=0;

//		ix=rand()%nx;
//		iy=rand()%ny;
//		iz=rand()%nz;
//configuration1		
		ix=3*nx/8;
		iy=ny/2;
		iz=nz/2;
		jx=ix+nx/4;
		jy=iy;
		jz=iz;
//configuration2		
//		ix=3*nx/8;
//		iy=ny/2;
//		iz=nz/2;
//		jx=ix+nx/4;
//		jy=iy-ny/8;
//		jz=iz;
//configuration3
//		ix=nx/2;
//		iy=ny/2+ny/8;
//		iz=nz/2-nz/8;
//		jx=ix+nx/8;
//		jy=iy-ny/4;
//		jz=iz+nz/4;

//		ix=iy=iz=nx/2;

//		for(i=0;i<nx;i++)
//			for(j=0;j<ny;j++)
//				for(k=0;k<nz;k++)
//				{
					indx2=index2(i,j,k);

					tp1=(float)(i-ix);
					if(fabs(tp1)>nx/2)
						tp1=tp1>0?(tp1-nx):(tp1+nx);
					tp2=(float)(j-iy);
					if(fabs(tp2)>ny/2)
						tp2=tp2>0?(tp2-ny):(tp2+ny);
					tp3=(float)(k-iz);
					if(fabs(tp3)>nz/2)
						tp3=tp3>0?(tp3-nz):(tp3+nz);
					mp1=(float)(i-jx);
					if(fabs(mp1)>nx/2)
						mp1=mp1>0?(mp1-nx):(mp1+nx);
					mp2=(float)(j-jy);
					if(fabs(mp2)>ny/2)
						mp2=mp2>0?(mp2-ny):(mp2+ny);
					mp3=(float)(k-jz);
					if(fabs(mp3)>nz/2)
						mp3=mp3>0?(mp3-nz):(mp3+nz);
					if(tp1*tp1+tp2*tp2+tp3*tp3<=a[1]*a[1]/4.0 &&fabs(tp1*(0.9045)+tp2*(0.3015)+tp3*(-0.3015))<=a[2]/2.0)
//				            if(fabs(tp1*(0.9045)+tp2*(0.3015)+tp3*(-0.3015))<=a[2]/2.0)
							//if(fabs(tp2)<=a[2]/2.0)
							{
								eta[0][indx2]=1;
								eta[6][indx2]=0.54;//X_Ni
								eta[7][indx2]=0.28;//X_Hf
							}
//					continue;
						//place H phase at (1/4.1/4,1/4)
					//else if(mp1*mp1+mp2*mp2+mp3*mp3<=a[1]*a[1]/4.0 &&fabs(mp1*(0.9045)+mp2*(0.3015)+mp3*(0.3015))<=a[2]/2.0 )//variant2
					else if(mp1*mp1+mp2*mp2+mp3*mp3<=a[1]*a[1]/4.0 &&fabs(mp1*(0.9045)+mp2*(0.3015)+mp3*(-0.3015))<=a[2]/2.0 )//variant1
	//				else if(mp1*mp1+mp2*mp2+mp3*mp3<=a[1]*a[1]/4.0 &&fabs(mp1*(-0.3015)+mp2*(0.9045)+mp3*(-0.3015))<=a[2]/2.0 )//variant3
			//		else if(mp1*mp1+mp2*mp2+mp3*mp3<=a[1]*a[1]/4.0 &&fabs(mp1*(0.3015)+mp2*(-0.3015)+mp3*(0.9045)))<=a[2]/2.0 )//variant4
//						 if(fabs(mp1*(0.9045)+mp2*(0.3015)+mp3*(0.3015))<=a[2]/2.0)
							{
								eta[0][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
//					continue;
//					switch(iv)
//					{
//					case 0:
					//	if(tp1*tp1/a[2]/a[2]+tp2*tp2/a[1]/a[1]+tp3*tp3/a[1]/a[1]<=1.0/4.0)
//							if(fabs(tp1*(0.9045)+tp2*(0.3015)+tp3*(-0.3015))<=a[2]/2.0)
//							//if(fabs(tp2)<=a[2]/2.0)
//							{
//								eta[0][indx2]=1;
//								eta[6][indx2]=0.54;//X_Ni
//								eta[7][indx2]=0.28;//X_Hf
//							}
//						else if(fabs(mp1*(0.9045)+mp2*(0.3015)+mp3*(0.3015))<=a[2]/2.0)
//							{
//								eta[1][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						{
					//	if(tp1*tp1/a[2]/a[2]+tp2*tp2/a[1]/a[1]+tp3*tp3/a[1]/a[1]<=1.0/4.0)
//							if(fabs(tp1*(0.9045)+tp2*(-0.3015)+tp3*(0.3015))<=a[2]/2.0)
							//if(fabs(tp2)<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						break;
//					case 1:
//						if(ih)
//						{
//							if(fabs(tp1*(0.9045)+tp2*(0.3015)+tp3*(0.3015))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						else
//						{
//							if(fabs(tp1*(0.9045)+tp2*(-0.3015)+tp3*(-0.3015))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						break;
//					case 2:
//						if(ih)
//						{
//
//							if(fabs(tp1*(0.3015)+tp2*(0.9045)+tp3*(-0.3015))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						else
//						{
//							if(fabs(tp1*(-0.3015)+tp2*(0.9045)+tp3*(0.3015))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						break;
//					case 3:
//						if(ih)
//						{
//							if(fabs(tp1*(0.3015)+tp2*(0.9045)+tp3*(0.3015))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						else
//						{
//							if(fabs(tp1*(-0.3015)+tp2*(0.9045)+tp3*(-0.3015))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						break;
//					case 4:
//						if(ih)
//						{
//							if(fabs(tp1*(0.3015)+tp2*(-0.3015)+tp3*(0.9045))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						else
//						{
//							if(fabs(tp1*(-0.3015)+tp2*(0.3015)+tp3*(0.9045))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						break;
//					case 5:
//						if(ih)
//						{
//							if(fabs(tp1*(0.3015)+tp2*(0.3015)+tp3*(0.9045))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						else
//						{
//							if(fabs(tp1*(-0.3015)+tp2*(-0.3015)+tp3*(0.9045))<=a[2]/2.0)
//							{
//								eta[iv][indx2]=1;
//								eta[6][indx2]=0.54;
//								eta[7][indx2]=0.28;
//							}
//						}
//						break;
//
//					default:
//						break;
						
//					}//end of switch iv
//				   }//end of space choosing

//				   }//end of particle number

				/*
				for(g1=0;g1<ng;g1++)
				{
					if(g1==gs[indx2])
					{
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);

							switch(v1)
							{
							case 0:
								if((i-nx/4)*(i-nx/4)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
									(eta[indx1])[indx2]=1;
								else
									(eta[indx1])[indx2]=0;
								break;
							case 1:
								if((i-nx*3/4)*(i-nx*3/4)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
									(eta[indx1])[indx2]=1;
								else
									(eta[indx1])[indx2]=0;
								break;
							case 4:
								if((i-nx/4)*(i-nx/4)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0] 
									|| (i-nx*3/4)*(i-nx*3/4)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
									(eta[indx1])[indx2]=0.375;
								else
									(eta[indx1])[indx2]=0.225;
								break;

							case 5:
								if((i-nx/4)*(i-nx/4)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0] 
									|| (i-nx*3/4)*(i-nx*3/4)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
									(eta[indx1])[indx2]=0.167;
								else
									(eta[indx1])[indx2]=0.233;
								break;
								
							default:
								eta[indx1][indx2]=0;
							}

				}
					}
					else
					{
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
							(eta[indx1])[indx2]=0;
						}
					}
				}
*/				
			}
		
	return 0;
}
// Two single phase particle 
int GInitConfig::Init_DOUBLE(float *eta[],float a[],int *gs)
{
	int g1,v1,i,j,k;
	long int indx1,indx2;

	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1==gs[indx2])
					{
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);

							switch(v1)
							{
							case 0:
								if(i<nx/2)
									(eta[indx1])[indx2]=1;
								else
									(eta[indx1])[indx2]=0;
								break;

							case 4:
								if(i<nx/2)
									(eta[indx1])[indx2]=0.534;
								else
									(eta[indx1])[indx2]=0.505;
								break;

							case 5:
								if(i<nx/2)
									(eta[indx1])[indx2]=0.28;
								else
									(eta[indx1])[indx2]=0.19;
								break;
								
							default:
								eta[indx1][indx2]=0;
							}

				}
					}
					else
					{
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
							(eta[indx1])[indx2]=0;
						}
					}
				}
			}
		
	return 0;
}

int GInitConfig::Init_TWIN_13(float *eta[],float a[],int *gs)
{
	int g1,i,j,k;
	long int indx1_1,indx1_3,indx2;

	if(nv<3)
	{
		cout<<"Invalid InitConfig option!"<<endl;
		1;
	}
	Init_HOMO(eta);
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
					indx1_1=index1(g1,0);
					indx1_3=index1(g1,2);
//	habit=(0.6969,0.6646,-0.2697) twin=(0,1,-1)
					if(i>=nx*3/16 && i<nx*13/16 && j>=ny*3/16.0 && j<ny*13/16.0 && k>=nz*3/16.0 && k<nz*13/16.0)
					{
					   	if(fabs((float)(i-nx/2)*(0.6969)+(float)(j-ny/2)*(0.6646)+(float)(k-nz/2)*(-0.2697))<=a[1])
						{
							if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*3.135))
							{
								(eta[indx1_1])[indx2]=0;
								(eta[indx1_3])[indx2]=1;
								if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*2.865))
								{
									(eta[indx1_1])[indx2]=1;
									(eta[indx1_3])[indx2]=0;
									if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*2.135))
									{
										(eta[indx1_1])[indx2]=0;
										(eta[indx1_3])[indx2]=1;
										if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*1.865))
										{
											(eta[indx1_1])[indx2]=1;
											(eta[indx1_3])[indx2]=0;
											if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*1.135))
											{
												(eta[indx1_1])[indx2]=0;
												(eta[indx1_3])[indx2]=1;
												if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*0.865))
												{
													(eta[indx1_1])[indx2]=1;
													(eta[indx1_3])[indx2]=0;
					 								if(fabs((float)(i-nx/2)*0+(float)(j-ny/2)*1+(float)(k-nz/2)*(-1))<=(a[2]*0.135))
													{
														(eta[indx1_1])[indx2]=0;
														(eta[indx1_3])[indx2]=1;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
	return 0;
}

// For cubic->tetragonal system
int GInitConfig::Init_TWIN_12(float *eta[],float a[],int *gs)
{
	int g1,i,j,k;
	long int indx1_1,indx1_2,indx1_3,indx1_4,indx2;

	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
					indx1_1=index1(g1,0);
					indx1_2=index1(g1,3);
					indx1_3=index1(g1,9);
					indx1_4=index1(g1,10);

//	twin=(110), no habit prefer
/*		
		if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
					{
						if((float)(i-nx/2)*1+(float)(j-ny/2)*0+(float)(k-nz/2)*0 <=(a[2]*1.25))
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							if(fabs((float)(i-nx/2)*1+(float)(j-ny/2)*1+(float)(k-nz/2)*0)<=(a[2]*0.75))
							{
								eta[indx1_1][indx2]=0;
								eta[indx1_2][indx2]=1;
								if(fabs((float)(i-nx/2)*1+(float)(j-ny/2)*1+(float)(k-nz/2)*0)<=(a[2]*0.25))
								{
									eta[indx1_1][indx2]=1;
									eta[indx1_2][indx2]=0;
								}
							}
						}
					}
*/
//twin structure of variant 1,4,10,11
   if(fabs(j-ny/2)<=a[0]/2)
       {
       //if(fabs(k-nz/2+a[0]/2)<=a[0]/2 && fabs(i+k+a[0]/2)<=a[0]/2)
       if(fabs(k-nz/2+a[0]/2)<=a[0]/2 && fabs(k-i+a[0]/2)<=a[0]/2)
       {
       eta[indx1_1][indx2]=1;
       eta[indx1_2][indx2]=0;
       eta[indx1_3][indx2]=0;
       eta[indx1_4][indx2]=0;
       }
       if(fabs(k-nz/2+a[0]/2)<=a[0]/2 && fabs(k-i-a[0]/2)<=a[0]/2)
       {
       eta[indx1_2][indx2]=1; 
       eta[indx1_1][indx2]=0;
       eta[indx1_3][indx2]=0;
       eta[indx1_4][indx2]=0;
       }
       if(fabs(k-nz/2-a[0]/2)<=a[0]/2 && fabs(i+k-nz-a[0]/2)<=a[0]/2)
       {
       eta[indx1_3][indx2]=0;
       eta[indx1_2][indx2]=0;
       eta[indx1_1][indx2]=0;
       eta[indx1_4][indx2]=1;
       }
       if(fabs(k-nz/2-a[0]/2)<=a[0]/2 && fabs(i+k-nz+a[0]/2)<=a[0]/2)
      {
      eta[indx1_4][indx2]=0;
       eta[indx1_2][indx2]=0;
       eta[indx1_3][indx2]=1;
       eta[indx1_1][indx2]=0;
      }
      }
				}
			}
	return 0;
}

//	For cubic to monoclinic TiNiPt prep
int GInitConfig::Init_TWIN_159(float *eta[],float a[],int *gs)
{
	int g1,i,j,k;
	long int indx1_1,indx1_2,indx1_3,indx2;
	float temp;

	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
					indx1_1=index1(g1,0);
					indx1_2=index1(g1,4);
					indx1_3=index1(g1,8);
//	twin=(100), no habit prefer
					if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<=a[0]*a[0])
					{
						temp=(float)(i-nx/2)*1+(float)(j-ny/2)*0+(float)(k-nz/2)*0;
						if(temp<=a[2]*2.5 && temp>a[2]*1.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*1.5 && temp>a[2]*0.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*0.5 && temp>a[2]*(-0.5))
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*(-0.5) && temp>a[2]*(-1.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*(-1.5) && temp>a[2]*(-2.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}
					}
				}
			}
	return 0;
}

// Load averaged 1 and redefine it as laminate structure piled on (111)
int GInitConfig::Init_LAM_123(float *eta[],char eta_file[],float a[],int *gs)
{
	int g1,i,j,k;
	long int indx1_1,indx1_2,indx1_3,indx2;
	float temp;
	float threshold=0.001;
	
	Init_LOAD(eta,eta_file);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
					indx1_1=index1(g1,0);
					indx1_2=index1(g1,1);
					indx1_3=index1(g1,2);

					if(eta[0][indx2]<threshold)
					{
						eta[indx1_1][indx2]=0;
						eta[indx1_2][indx2]=0;
						eta[indx1_3][indx2]=0;
						continue;
					}
					else
					{
//	twin=(111)
						temp=(float)(i-nx/2)*0.5733+(float)(j-ny/2)*0.5733+(float)(k-nz/2)*0.5733;
						if(temp<=a[2]*10.5 && temp>a[2]*9.5)
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*9.5 && temp>a[2]*8.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*8.5 && temp>a[2]*7.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}

						if(temp<=a[2]*7.5 && temp>a[2]*6.5)
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*6.5 && temp>a[2]*5.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*5.5 && temp>a[2]*4.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}

						if(temp<=a[2]*4.5 && temp>a[2]*3.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*3.5 && temp>a[2]*2.5)
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*2.5 && temp>a[2]*1.5)
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}

						if(temp<=a[2]*1.5 && temp>a[2]*0.5)
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*0.5 && temp>a[2]*(-0.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*(-0.5) && temp>a[2]*(-1.5))
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}

						if(temp<=a[2]*(-1.5) && temp>a[2]*(-2.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*(-2.5) && temp>a[2]*(-3.5))
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*(-3.5) && temp>a[2]*(-4.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}

						if(temp<=a[2]*(-4.5) && temp>a[2]*(-5.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*(-5.5) && temp>a[2]*(-6.5))
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*(-6.5) && temp>a[2]*(-7.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}

						if(temp<=a[2]*(-7.5) && temp>a[2]*(-8.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=1;
							eta[indx1_3][indx2]=0;
						}
						if(temp<=a[2]*(-8.5) && temp>a[2]*(-9.5))
						{
							eta[indx1_1][indx2]=0;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=1;
						}
						if(temp<=a[2]*(-9.5) && temp>a[2]*(-10.5))
						{
							eta[indx1_1][indx2]=1;
							eta[indx1_2][indx2]=0;
							eta[indx1_3][indx2]=0;
						}
					}
				}
			}

	return 0;
}

int GInitConfig::Init_USER(float *eta[],float a[],int *gs,float *conc1,float *conc2)
{
	Init_HOMO(eta);
	char filename[500];
	
//	sprintf(filename,"Mplate_v1.vtk");
//	sprintf(filename,"%s",input_eta);
	//init_from_file(*eta,filename);

//	sprintf(filename,"input_v1.vtk");
//	sprintf(filename,"V1_input.vtk");
	//sprintf(filename,"M_p3_2_v1_49.vtk");
	//sprintf(filename,"pn0_199.vtk");
	//sprintf(filename,"M_pure_7_v1_re59.vtk");
	//sprintf(filename,"Big2_r_v1_99.vtk");
//	sprintf(filename,"ther_T02_r_v1_49.vtk");
//	sprintf(filename,"test_2d_v1_49.vtk");
	//sprintf(filename,"ther_str_v1_39.vtk");
	//sprintf(filename,"tetra_v1_3_27.vtk");
//	sprintf(filename,"test_stress_T04_v1_49.vtk");

//	sprintf(filename,"Small2_r_v1_99.vtk");
//	sprintf(filename,"ther_1_v1_49.vtk");
	//sprintf(filename,"test3_v1_99.vtk");
//	sprintf(filename,"Stress_par_1_v1_27.vtk");
//	sprintf(filename,"stress_test3_v1_19.vtk");

//	sprintf(filename,"stress_c1_v1_19.vtk");
	//sprintf(filename,"stress_e4_T190_v1_32.vtk");
/*	sprintf(filename,"stress_e4c_T200_v1_19.vtk");
//	sprintf(filename,"stress_test_T180_v1_19.vtk");
//	sprintf(filename,"eta_v1_rand_18.vtk");
	init_from_file(&eta[0][0],filename);
//	sprintf(filename,"T160_e4r_v2_19.vtk");
//	sprintf(filename,"stress_c1_v2_19.vtk");
//	sprintf(filename,"stress_t_T180_v2_19.vtk");
//	sprintf(filename,"stress_e4_T190_v2_32.vtk");
//	sprintf(filename,"stress_test_T180_v2_19.vtk");
	sprintf(filename,"stress_e4c_T200_v2_19.vtk");
	init_from_file(&eta[1][0],filename);
//	sprintf(filename,"T160_e4r_v3_19.vtk");
//	sprintf(filename,"stress_c1_v3_19.vtk");
//	sprintf(filename,"stress_t_T180_v3_19.vtk");
//	sprintf(filename,"stress_e4_T190_v3_32.vtk");
//	sprintf(filename,"stress_test_T180_v3_19.vtk");
	sprintf(filename,"stress_e4c_T200_v3_19.vtk");
	init_from_file(&eta[2][0],filename);
//	sprintf(filename,"T160_e4r_v4_19.vtk");
//	sprintf(filename,"stress_c1_v4_19.vtk");
//	sprintf(filename,"stress_t_T180_v4_19.vtk");
//	sprintf(filename,"stress_e4_T190_v4_32.vtk");
//	sprintf(filename,"stress_test_T180_v4_19.vtk");
	sprintf(filename,"stress_e4c_T200_v4_19.vtk");
	init_from_file(&eta[3][0],filename);
*/
/*	sprintf(filename,"T160_c1_v5_08.vtk");
	init_from_file(&eta[4][0],filename);
	sprintf(filename,"T160_c1_v6_08.vtk");
	init_from_file(&eta[5][0],filename);
	sprintf(filename,"T160_c1_v7_08.vtk");
	init_from_file(&eta[6][0],filename);
	sprintf(filename,"T160_c1_v8_08.vtk");
	init_from_file(&eta[7][0],filename);
	sprintf(filename,"T160_c1_v9_08.vtk");
	init_from_file(&eta[8][0],filename);
	sprintf(filename,"T160_c1_v10_08.vtk");
	init_from_file(&eta[9][0],filename);
	sprintf(filename,"T160_c1_v11_08.vtk");
	init_from_file(&eta[10][0],filename);
	sprintf(filename,"T160_c1_v12_08.vtk");
	init_from_file(&eta[11][0],filename);
*/
//sprintf(filename,"Mplate_v2.vtk");
	//sprintf(filename,"M_p3_2_v2_49.vtk");
//	sprintf(filename,"V2_input.vtk");
//	sprintf(filename,"ther_1_v2_49.vtk");
//	sprintf(filename,"Small2_v2_99.vtk");
//	sprintf(filename,"test_2d_v2_49.vtk");
//	sprintf(filename,"ther_T02_r_v2_49.vtk");
//	sprintf(filename,"test_stress_T04_v2_49.vtk");
//	sprintf(filename,"Mart_plat4_v2_13.vtk");
//	sprintf(filename,"Big2_r_v2_99.vtk");
	//sprintf(filename,"Stress_par_1_v2_27.vtk");
	//sprintf(filename,"pn1_199.vtk");
	//sprintf(filename,"input_v2.vtk");
	//sprintf(filename,"ther_str_v2_39.vtk");
	//sprintf(filename,"tetra_v2_3_27.vtk");
//	sprintf(filename,"stress_test3_v2_19.vtk");
//	sprintf(filename,"stress_d1_v2_19.vtk");
	//sprintf(filename,"tetra_v2_15.vtk");
//	init_from_file(&eta[1][0],filename);
//	sprintf(filename,"Mplate_v3.vtk");
	//sprintf(filename,"M_p3_2_v3_49.vtk");
//	sprintf(filename,"V3_input.vtk");
//	sprintf(filename,"ther_1_v3_49.vtk");
//	sprintf(filename,"Small2_v3_99.vtk");
//	sprintf(filename,"test_2d_v3_49.vtk");
//	sprintf(filename,"Mart_plat4_v3_13.vtk");
//	sprintf(filename,"Big2_r_v3_99.vtk");
	//sprintf(filename,"Stress_par_1_v3_27.vtk");
//	sprintf(filename,"pn2_199.vtk");
	//sprintf(filename,"input_v3.vtk");
//	sprintf(filename,"ther_T02_r_v3_49.vtk");
	//sprintf(filename,"ther_str_v3_39.vtk");
//	sprintf(filename,"stress_d1_v3_19.vtk");
	//sprintf(filename,"stress_test3_v3_19.vtk");
	//sprintf(filename,"tetra_v3_15.vtk");
//	init_from_file(&eta[2][0],filename);
//	sprintf(filename,"Mplate_v4.vtk");
//	sprintf(filename,"stress_d1_v4_19.vtk");
	//sprintf(filename,"stress_test3_v4_19.vtk");
	//sprintf(filename,"tetra_v3_15.vtk");
//	init_from_file(&eta[3][0],filename);
	sprintf(filename,"eta_v1_rand3_12.vtk");
	init_from_file(&eta[4][0],filename);

	sprintf(filename,"eta_v2_rand3_12.vtk");
	init_from_file(&eta[5][0],filename);
	sprintf(filename,"eta_v3_rand3_12.vtk");
	init_from_file(&eta[6][0],filename);
	sprintf(filename,"eta_v4_rand3_12.vtk");
	init_from_file(&eta[7][0],filename);
	sprintf(filename,"eta_v5_rand3_12.vtk");
	init_from_file(&eta[8][0],filename);
	sprintf(filename,"eta_v6_rand3_12.vtk");
	init_from_file(&eta[9][0],filename);

//sprintf(filename,"M_p3_2_v4_49.vtk");
//	sprintf(filename,"V4_input.vtk");
//	sprintf(filename,"Small2_v4_99.vtk");
//	sprintf(filename,"ther_1_v4_49.vtk");
//	sprintf(filename,"ther_T02_v4_49.vtk");
//	sprintf(filename,"Mart_plat4_v4_13.vtk");
//	sprintf(filename,"Big2_r_v4_99.vtk");
//	sprintf(filename,"Stress_par_1_v4_27.vtk");
//	sprintf(filename,"pn3_199.vtk");
	//sprintf(filename,"input_v4.vtk");/
//	sprintf(filename,"ther_T02_r_v4_49.vtk");
//	sprintf(filename,"ther_str_v4_39.vtk");
//	init_from_file(&eta[3][0],filename);
	
//	sprintf(filename,"Hphase_p3_v1_02.vtk");
//	sprintf(filename,"%s",input_eta);
	//init_from_file(*eta,filename);
//	sprintf(filename,"Plate_v1.vtk");
	sprintf(filename,"Ni_rand3_12.vtk");
//	sprintf(filename,"Ni_rand_12.vtk");
//	sprintf(filename,"Hphase_large2_Ni_10.vtk");
//	//sprintf(filename,"Hphase_1d_Ni_49.vtk");
//	init_from_file(&eta[6][0],filename);
	init_from_file(&conc1[0],filename);
	sprintf(filename,"Hf_rand3_12.vtk");
//	sprintf(filename,"Hf_diss_18.vtk");
//	sprintf(filename,"Hf_conc2_05.vtk");
//	sprintf(filename,"Hphase_large2_Hf_10.vtk");
//	init_from_file(&eta[7][0],filename);
	init_from_file(&conc2[0],filename);

//	sprintf(filename,"Part_V1.vtk");
//	init_from_file(&eta[4][0],filename);
//	sprintf(filename,"Hphase_p3_v2_02.vtk");
//	sprintf(filename,"Part_V2.vtk");
//	sprintf(filename,"Plate_v2.vtk");
//	sprintf(filename,"Hphase_p_v2_09.vtk");
//	sprintf(filename,"input_v2.vtk");
//	init_from_file(&eta[5][0],filename);
//	sprintf(filename,"Hphase_p_v3_09.vtk");
//	sprintf(filename,"Part_V3.vtk");
//	sprintf(filename,"Hphase_p3_v3_02.vtk");
//	sprintf(filename,"Plate_v3.vtk");
	//sprintf(filename,"hph_v3_01.vtk");
//	sprintf(filename,"input_v3.vtk");
//	init_from_file(&eta[6][0],filename);
//	sprintf(filename,"Hphase_p_v4_09.vtk");
//	sprintf(filename,"Plate_v4.vtk");
//	sprintf(filename,"Part_V4.vtk");
	//sprintf(filename,"Hphase_p3_v4_02.vtk");
//	sprintf(filename,"input_v4.vtk");
//	init_from_file(&eta[7][0],filename);
//	sprintf(filename,"Hphase_p_v5_09.vtk");
//	sprintf(filename,"Part_V5.vtk");
//	sprintf(filename,"Hphase_p3_v5_02.vtk");
//	sprintf(filename,"Plate_v5.vtk");
//	sprintf(filename,"input_v5.vtk");
//	init_from_file(&eta[8][0],filename);
//	sprintf(filename,"Hphase_p_v6_09.vtk");
//	sprintf(filename,"Plate_v6.vtk");
//	sprintf(filename,"Part_V6.vtk");
	//sprintf(filename,"Hphase_p3_v6_02.vtk");
//	sprintf(filename,"input_v6.vtk");
//	init_from_file(&eta[9][0],filename);

/*
	sprintf(filename,"PlatV5_10.vtk");
//	sprintf(filename,"input_v5.vtk");
	init_from_file(&eta[8][0],filename);
	sprintf(filename,"PlatV6_10.vtk");
//	sprintf(filename,"input_v6.vtk");
	init_from_file(&eta[9][0],filename);
	*/
/*	
	sprintf(filename,"Plate_v1.vtk");
//        sprintf(filename,"hph_v1_01.vtk");
//	sprintf(filename,"Hphase_p_v1_09.vtk");
	init_from_file(&eta[4][0],filename);
	sprintf(filename,"Plate_v2.vtk");
//	sprintf(filename,"Hphase_p_v2_09.vtk");
//	sprintf(filename,"input_v2.vtk");
//	sprintf(filename,"hph_v2_01.vtk");
	init_from_file(&eta[5][0],filename);
//	sprintf(filename,"hph_v3_01.vtk");
//	sprintf(filename,"Hphase_p_v3_09.vtk");
	sprintf(filename,"Plate_v3.vtk");
	//sprintf(filename,"hph_v3_01.vtk");
//	sprintf(filename,"input_v3.vtk");
	init_from_file(&eta[6][0],filename);
//	sprintf(filename,"Hphase_p_v4_09.vtk");
	sprintf(filename,"Plate_v4.vtk");
//	sprintf(filename,"hph_v4_01.vtk");
//	sprintf(filename,"input_v4.vtk");
	init_from_file(&eta[7][0],filename);
//	sprintf(filename,"hph_v5_01.vtk");
//	sprintf(filename,"Hphase_p_v5_09.vtk");
	sprintf(filename,"Plate_v5.vtk");
//	sprintf(filename,"input_v5.vtk");
	init_from_file(&eta[8][0],filename);
//	sprintf(filename,"Hphase_p_v6_09.vtk");
	sprintf(filename,"Plate_v6.vtk");
//	sprintf(filename,"hph_v6_01.vtk");
	init_from_file(&eta[9][0],filename);
*/
return 0;
}

int GInitConfig::Init_LOAD(float *eta[],char eta_file[])
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	char temp[80]={0};

	if(eta_file==NULL)
	{
		Init_HOMO(eta);
		return 1;
	}

	ifstream fin(eta_file,ios::in);
	if(!fin)
	{
		cerr<<"Eta input can not be found!";
		1;
	}
	for(i=0;i<30;i++) // ignore VTK header
		fin>>temp;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						fin>>(eta[indx1])[indx2];
					}
		}
	fin.close();

	return 0;
}

int GInitConfig::Init_LOAD_3264C(float *eta[],char eta_file[])
{
	int g1,v1,i,j,k;
	long int indx1,indx2;
	char temp[80]={0};
	float tp1;
	int no1;

	if(eta_file==NULL)
	{
		Init_HOMO(eta);
		return 1;
	}

	ifstream fin(eta_file,ios::in);
	if(!fin)
	{
		cerr<<"Eta input can not be found!";
		1;
	}
	for(i=0;i<30;i++) // ignore VTK header
		fin>>temp;
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			tp1=0;
			no1=0;
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						if(i>=nx/4 && i<nx*3/4 && j>=ny/4 && j<ny*3/4 && k>=nz/4 && k<nz*3/4)
						{
							fin>>eta[indx1][indx2];
							if(i==nx/4 || i==nx*3/4-1 || j==ny/4 || j==j*3/4-1 || k==nz/4 || k==nz*3/4-1)
								switch(v1)
								{
								case 0:
								case 1:
								case 2:
								case 3:
									break;
								case 4:
								case 5:
									tp1+=eta[indx1][indx2];
									no1++;
									break;
								}
						}
						else
						{
							switch(v1)
							{
							case 0:
							case 1:
							case 2:
							case 3:
								eta[indx1][indx2]=0;
								break;
							case 4:
								eta[indx1][indx2]=0.285;
								break;
							case 5:
								eta[indx1][indx2]=0.207;
								break;
							}
						}
					}
			
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						if(i>=nx/4 && i<nx*3/4 && j>=ny/4 && j<ny*3/4 && k>=nz/4 && k<nz*3/4)
							continue;
						switch(v1)
						{
						case 0:
						case 1:
						case 2:
						case 3:
							eta[indx1][indx2]=0;
							break;
						case 4:
						case 5:
							eta[indx1][indx2]=tp1/no1;
							break;
						}
					}
			

		}
	fin.close();

	return 0;
}

int GInitConfig::Init_Hphase(float *eta[],float a[],int *gs) // a[0]=number of particles; a[1]=particle length; a[2]=width;
{
	int i,j,k;
	int ip,iv,ix,iy,iz,ih,ix_old,iy_old,iz_old;
	long int indx2;
	float tp1,tp2,tp3,t1,t2,t3,m1,m2,m3;

	srand(120);

	//srand(0);

	Init_HOMO(eta);

		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)		
		{
					indx2=index2(i,j,k);
		if(k<=nz/4)
		{
							eta[0][indx2]=1;
							eta[6][indx2]=0.53;//x_Ni
						eta[7][indx2]=0.26;//x-Hf
}
}
/*
else if((k>=nz/2)&&(k<=3*nz/4))
{
							eta[0][indx2]=0;
							eta[6][indx2]=0.498;//x_Ni
						eta[7][indx2]=0.19;//x-Hf
}
}
*/
//				}
//	For 2D strain field of H-phase
/*
	for (i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				if(fabs((i-nx/2)*1.414+(j-ny/2)*3)<a[1]
					&& fabs((i-nx/2)*(-1.414)+(j-ny/2)*3)<a[1])
					{
						eta[6][indx2]=1;
						eta[4][indx2]=0.375;
						eta[5][indx2]=0.167;
					}
			}
*/
//	For short time particle with {100} habit
//
/*
	for(ip=0;ip<(int)a[0];ip++)
	{
	//	iv=rand()%6;
		iv=0;

//		ix=rand()%nx;
//		iy=rand()%ny;
//		iz=rand()%nz;
		

                ix=nx/2;
		iy=ny/2;
		iz=nz/2;
		
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);

					tp1=fabs((float)(i-ix));
					tp1=tp1<nx/2?tp1:nx-tp1;
					tp2=fabs((float)(j-iy));
					tp2=tp1<ny/2?tp2:ny-tp2;
					tp3=fabs((float)(k-iz));
					tp3=tp3<nz/2?tp3:nz-tp3;

					switch(iv)
					{
					case 0:
					case 1:
						if(tp1*tp1/a[2]/a[2]+tp2*tp2/a[1]/a[1]+tp3*tp3/a[1]/a[1]<=1.0/4.0)
						{
							eta[iv][indx2]=1;
							eta[6][indx2]=0.59;//x_Ni
						eta[7][indx2]=0.29;//x-Hf
						}
						break;
					case 2:
					case 3:
						if(tp1*tp1/a[1]/a[1]+tp2*tp2/a[2]/a[2]+tp3*tp3/a[1]/a[1]<=1.0/4.0)
						{
							eta[iv][indx2]=1;
							eta[6][indx2]=0.53;//x_Ni
							eta[7][indx2]=0.268;//x_Hf
						}
						break;
					case 4:
					case 5:
						if(tp1*tp1/a[1]/a[1]+tp2*tp2/a[1]/a[1]+tp3*tp3/a[2]/a[2]<=1.0/4.0)
						{
							eta[iv][indx2]=1;
							eta[6][indx2]=0.53;//x_Ni
							eta[7][indx2]=0.268;//x_Hf
						}
						break;
					default:
						break;
					}
				}

}
*/
//	For long time particle with {311} habit
/*
	cout<<"Generation of H-phase particles with {311} habit plane."<<endl;
//	ix_old=0;
//	iy_old=0;
//	iz_old=0;
	for(ip=0;ip<(int)a[0];ip++)
	{
//       		iv=rand()%6;
		iv=0;
	//iv=0;
//		ih=rand()%2;
//		ih=0;

//		ix=rand()%nx;
//		iy=rand()%ny;
//		iz=rand()%nz;
					//if((ix-ix_old)*(ix-ix_old)+(iy-iy_old)*(iy-iy_old)+(iz-iz_old)*(iz-iz_old)<=a[1]*a[1])
					//	continue;

		ix=iy=iz=nx/2;

		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);

					tp1=(float)(i-ix);
					if(fabs(tp1)>nx/2)
						tp1=tp1>0?(tp1-nx):(tp1+nx);
					tp2=(float)(j-iy);
					if(fabs(tp2)>ny/2)
						tp2=tp2>0?(tp2-ny):(tp2+ny);
					tp3=(float)(k-iz);
					if(fabs(tp3)>nz/2)
						tp3=tp3>0?(tp3-nz):(tp3+nz);
						//place H phase at (1/4.1/4,1/4)
					if(tp1*tp1+tp2*tp2+tp3*tp3>a[1]*a[1]/4.0)
					continue;
					switch(iv)
					{
					case 0:
						if(ih)
						{
					//	if(tp1*tp1/a[2]/a[2]+tp2*tp2/a[1]/a[1]+tp3*tp3/a[1]/a[1]<=1.0/4.0)
							if(fabs(tp1*(0.9045)+tp2*(0.3015)+tp3*(-0.3015))<=a[2]/2.0)
							//if(fabs(tp2)<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;//X_Ni
								eta[7][indx2]=0.28;//X_Hf
							}
						}
						else
						{
					//	if(tp1*tp1/a[2]/a[2]+tp2*tp2/a[1]/a[1]+tp3*tp3/a[1]/a[1]<=1.0/4.0)
							if(fabs(tp1*(0.9045)+tp2*(-0.3015)+tp3*(0.3015))<=a[2]/2.0)
							//if(fabs(tp2)<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						break;
					case 1:
						if(ih)
						{
							if(fabs(tp1*(0.9045)+tp2*(0.3015)+tp3*(0.3015))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						else
						{
							if(fabs(tp1*(0.9045)+tp2*(-0.3015)+tp3*(-0.3015))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						break;
					case 2:
						if(ih)
						{

							if(fabs(tp1*(0.3015)+tp2*(0.9045)+tp3*(-0.3015))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						else
						{
							if(fabs(tp1*(-0.3015)+tp2*(0.9045)+tp3*(0.3015))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						break;
					case 3:
						if(ih)
						{
							if(fabs(tp1*(0.3015)+tp2*(0.9045)+tp3*(0.3015))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						else
						{
							if(fabs(tp1*(-0.3015)+tp2*(0.9045)+tp3*(-0.3015))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						break;
					case 4:
						if(ih)
						{
							if(fabs(tp1*(0.3015)+tp2*(-0.3015)+tp3*(0.9045))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						else
						{
							if(fabs(tp1*(-0.3015)+tp2*(0.3015)+tp3*(0.9045))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						break;
					case 5:
						if(ih)
						{
							if(fabs(tp1*(0.3015)+tp2*(0.3015)+tp3*(0.9045))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						else
						{
							if(fabs(tp1*(-0.3015)+tp2*(-0.3015)+tp3*(0.9045))<=a[2]/2.0)
							{
								eta[iv][indx2]=1;
								eta[6][indx2]=0.54;
								eta[7][indx2]=0.28;
							}
						}
						break;

					default:
						break;
						
					}//end of switch iv
				   }//end of space choosing

				   }//end of particle number
*/
return 0;
}
 
void GInitConfig::Init_field(float *eta[],float *conc1,float *conc2)
{
	Init_HOMO(eta);
	char filename[500];
	
//	sprintf(filename,"Plate_v1.vtk");
//	sprintf(filename,"Part_V1.vtk");
	//sprintf(filename,"Hphase_p_v1_00.vtk");
//	sprintf(filename,"Hphase_two2_v1_05.vtk");
	sprintf(filename,"eta_v1_rand5_10.vtk");
//	sprintf(filename,"Hphase_1d_v1_49.vtk");
//	sprintf(filename,"Hphase_p4_v1_02.vtk");
//	sprintf(filename,"%s",input_eta);
	//init_from_file(*eta,filename);
//	init_from_file(&eta[12][0],filename);
	init_from_file(&eta[4][0],filename);
//        sprintf(filename,"Part_V2.vtk");
	//sprintf(filename,"Hphase_p_v2_00.vtk");
	//sprintf(filename,"Hphase_small2.vtk");
//	sprintf(filename,"Hphase_two2_v2_05.vtk");
//	sprintf(filename,"Hphase_large2_v2_10.vtk");
//	sprintf(filename,"Hphase_1d_v2_49.vtk");
//	sprintf(filename,"Plate_v2.vtk");
//	sprintf(filename,"Hphase_p_v2_09.vtk");
//	sprintf(filename,"input_v2.vtk");
	sprintf(filename,"eta_v2_rand5_10.vtk");
	init_from_file(&eta[5][0],filename);
//	sprintf(filename,"Hphase_two2_v3_05.vtk");
//	sprintf(filename,"Hphase_large2_v3_10.vtk");
	//sprintf(filename,"Hphase_p_v3_00.vtk");
//	sprintf(filename,"Hphase_1d_v3_49.vtk");
//	sprintf(filename,"Hphase_small3.vtk");
	//sprintf(filename,"Hphase_p4_v3_02.vtk");
//	sprintf(filename,"Part_V3.vtk");
//	sprintf(filename,"Plate_v3.vtk");
	//sprintf(filename,"hph_v3_01.vtk");
//	sprintf(filename,"input_v3.vtk");
	sprintf(filename,"eta_v3_rand5_10.vtk");
	init_from_file(&eta[6][0],filename);
//	sprintf(filename,"Hphase_two2_v4_05.vtk");
//	sprintf(filename,"Hphase_large2_v4_10.vtk");
//	sprintf(filename,"Hphase_1d_v4_49.vtk");
//	sprintf(filename,"Hphase_small4.vtk");
//      sprintf(filename,"Hphase_p_v4_00.vtk");
//	sprintf(filename,"Plate_v4.vtk");
//	sprintf(filename,"Part_V4.vtk");
//	sprintf(filename,"Hphase_p4_v4_02.vtk");
//	sprintf(filename,"input_v4.vtk");
	sprintf(filename,"eta_v4_rand5_10.vtk");
	init_from_file(&eta[7][0],filename);
//	sprintf(filename,"Hphase_two2_v5_05.vtk");
//	sprintf(filename,"Hphase_large2_v5_10.vtk");
//	sprintf(filename,"Hphase_1d_v5_49.vtk");
//	sprintf(filename,"Hphase_small5.vtk");
//	sprintf(filename,"Hphase_p_v5_00.vtk");
	//sprintf(filename,"Hphase_p4_v5_02.vtk");
//	sprintf(filename,"Part_V5.vtk");
//	sprintf(filename,"Plate_v5.vtk");
//	sprintf(filename,"input_v5.vtk");
	sprintf(filename,"eta_v5_rand5_10.vtk");
	init_from_file(&eta[8][0],filename);
	//sprintf(filename,"Hphase_small6.vtk");
//	sprintf(filename,"Hphase_two2_v2_05.vtk");
//	sprintf(filename,"Hphase_large2_v6_10.vtk");
	//sprintf(filename,"Hphase_1d_v6_49.vtk");
//	sprintf(filename,"Hphase_p_v6_00.vtk");
//	sprintf(filename,"Plate_v6.vtk");
	//sprintf(filename,"Hphase_p4_v6_02.vtk");
//	sprintf(filename,"Part_V6.vtk");
//	sprintf(filename,"input_v6.vtk");
	sprintf(filename,"eta_v6_rand5_10.vtk");
	init_from_file(&eta[9][0],filename);
//	sprintf(filename,"Ni_conc2_05.vtk");
//	sprintf(filename,"Ni_diss_18.vtk");
	sprintf(filename,"Ni_rand5_10.vtk");
//	sprintf(filename,"Hphase_large2_Ni_10.vtk");
//	//sprintf(filename,"Hphase_1d_Ni_49.vtk");
//	init_from_file(&eta[6][0],filename);
	init_from_file(&conc1[0],filename);
	sprintf(filename,"Hf_rand5_10.vtk");
//	sprintf(filename,"Hf_diss_18.vtk");
//	sprintf(filename,"Hf_conc2_05.vtk");
//	sprintf(filename,"Hphase_large2_Hf_10.vtk");
//	init_from_file(&eta[7][0],filename);
	init_from_file(&conc2[0],filename);
 
//	sprintf(filename,"eta_input.vtk");
//	init_from_file(&eta[12][0],filename);
}
int     GInitConfig::Init_particle(float *eta[],float a[],int *gs)
{        
        int i,j,k;	
	int ip,iv,ix,iy,iz,ih,ix_old,iy_old,iz_old;
	long int indx2;
	float tp1,tp2,tp3,t1,t2,t3,m1,m2,m3;
	float radius;
	srand(100);
	random_seed=130;
	random_ave=8.1;
	random_var=1;
	


	Init_HOMO(eta);
	cout<<"Generation of sphere particles with Gaussian distribution size."<<endl;
	for(ip=0;ip<(int)a[0];ip++)
	{
	radius=Gauss_rand();
	cout<<radius/2<<endl;
	iv=0;
		ix=rand()%nx;
		iy=rand()%ny;
		iz=rand()%nz;
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);

					tp1=(float)(i-ix);
					if(fabs(tp1)>nx/2)
						tp1=tp1>0?(tp1-nx):(tp1+nx);
					tp2=(float)(j-iy);
					if(fabs(tp2)>ny/2)
						tp2=tp2>0?(tp2-ny):(tp2+ny);
					tp3=(float)(k-iz);
					if(fabs(tp3)>nz/2)
						tp3=tp3>0?(tp3-nz):(tp3+nz);
						//place H phase at (1/4.1/4,1/4)
				//	if(tp1*tp1+tp2*tp2+tp3*tp3>radius*radius/4.0)
				//	continue;
						if(tp1*tp1/radius/radius+tp2*tp2/radius/radius+tp3*tp3/radius/radius<=1.0/4.0)
	{	eta[iv][indx2]=1;
		eta[4][indx2]=0.54;//X_Ni
		eta[5][indx2]=0.28;//X_Hf
}
        }	//end of space loop
     }//end of particle number loop
        return 0;
}


int	GInitConfig::Init_GRAIN_3(float *eta[])
{
	int i,j,k;
	long int indx2;

	cout<<"This InitConfig is designed for generate 3 grains."<<endl;
	if(ng!=1 || nv!=3)
	{
		cout<<"Invalid ng or nv!"<<endl;
		1;
	}
	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				if((i-0)*(i-0)+(j-ny/2)*(j-ny/2)+(k-nz/2)*(k-nz/2)<5*5)
					(eta[0])[indx2]=1;
				if((i-nx/2)*(i-nx/2)+(j-0)*(j-0)+(k-nz/2)*(k-nz/2)<5*5)
					(eta[1])[indx2]=1;
				if((i-nx/2)*(i-nx/2)+(j-ny/2)*(j-ny/2)+(k-0)*(k-0)<5*5)
					(eta[2])[indx2]=1;
			}

	return 0;
}

int GInitConfig::Init_GRAIN_4(float *eta[])
{
	int i,j,k;
	long int indx2;

	cout<<"This InitConfig is designed for generate 4 grains."<<endl;
	if(ng!=1 || nv!=4)
	{
		cout<<"Invalid ng or nv!"<<endl;
		1;
	}
	Init_HOMO(eta);

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				if((i-nx/4)*(i-nx/4)+(j-ny/4)*(j-ny/4)+(k-nz/4)*(k-nz/4)<5*5)
					(eta[0])[indx2]=1;
				if((i-nx/4)*(i-nx/4)+(j-ny*3/4)*(j-ny*3/4)+(k-nz*3/4)*(k-nz*3/4)<5*5)
					(eta[1])[indx2]=1;
				if((i-nx*3/4)*(i-nx*3/4)+(j-ny/4)*(j-ny/4)+(k-nz*3/4)*(k-nz*3/4)<5*5)
					(eta[2])[indx2]=1;
				if((i-nx*3/4)*(i-nx*3/4)+(j-ny*3/4)*(j-ny*3/4)+(k-nz/4)*(k-nz/4)<5*5)
					(eta[3])[indx2]=1;

			}

	return 0;
}

int GInitConfig::Init_gs(int *gs,char gs_file[])
{
	int i,j,k;
	long int indx2;

	if(gs_file==NULL)
	{
	/*
		for(i=0;i<nx/2;i++)
			for(j=0;j<ny/2;j++)
				for(k=0;k<nz/2;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=0;
				}//grain 1
		for(i=nx/2;i<nx;i++)
			for(j=0;j<ny/2;j++)
				for(k=0;k<nz/2;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=1;
				}//grain 2
		for(i=nx/2;i<nx;i++)
			for(j=ny/2;j<ny;j++)
				for(k=0;k<nz/2;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=2;
				}//grain 3
		for(i=nx/2;i<nx;i++)
			for(j=ny/2;j<ny;j++)
				for(k=nz/2;k<nz;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=3;
				}//grain 4
		for(i=0;i<nx/2;i++)
			for(j=ny/2;j<ny;j++)
				for(k=0;k<nz/2;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=4;
				}//grain 5
		for(i=0;i<nx/2;i++)
			for(j=ny/2;j<ny;j++)
				for(k=nz/2;k<nz;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=5;
				}//grain 6
		for(i=0;i<nx/2;i++)
			for(j=0;j<ny/2;j++)
				for(k=nz/2;k<nz;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=6;
				}//grain 7
		for(i=nx/2;i<nx;i++)
			for(j=0;j<ny/2;j++)
				for(k=nz/2;k<nz;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=7;
				}//grain 8
*/
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);
					gs[indx2]=0;
				}//grain 1
		return 1;
	}
        
	ifstream fin(gs_file,ios::in);

	if(!fin)
	{
		cerr<<"Grain shape input can not be found!";
		1;
	}

	for(k=0;k<nz;k++)
		for(j=0;j<ny;j++)
			for(i=0;i<nx;i++)
			{
				indx2=index2(i,j,k);
				fin>>gs[indx2];
			}
	fin.close();

	return 0;
	
}

void GInitConfig::init_from_file(float *tpr, char *filename)
{
	int x_i,x_j,x_k,xijk,i;
	int l;
	FILE *fp=NULL;
	char ignore[40];
	int fnx,fny,fnz;
fp=fopen(filename,"r");
	if(fp)
	{
		l=strlen(filename);
		if(strcmp(&filename[l-4],".dat")==0)
		{	for_xijk
			fscanf(fp,"%lf",&tpr[xijk]);
			efor_xijk
		}
		else if(strcmp(&filename[l-4],".vtk")==0)
		{	for(i=0;i<4;i++)
				fgets(ignore,40,fp);//ignoring first four lines	
			fscanf(fp,"DIMENSIONS %d %d %d\n",&fnx,&fny,&fnz);
			for(i=0;i<5;i++)
				fgets(ignore,40,fp);//ignoring 5 lines
			if( (nx!=fnx) | (ny!=fny) | (nz!=fnz) )
				{fprintf(stderr,"#Error:Input file size (%d,%d,%d) is different from the expected size (%d,%d,%d)\n",fnx,fny,fnz,nx,ny,nz);exit(0);}
for(x_k=0; x_k<nz; x_k++)
	for(x_j=0; x_j<ny; x_j++)
		for(x_i=0; x_i<nx; x_i++)
			//fscanf(fp,"%lf\n",&tpr[x_i*nx*ny+x_j*ny+x_k]);
			fscanf(fp,"%f",&tpr[x_i*ny*nz+x_j*nz+x_k]);
                         printf("Input from user set succeed.\n");
		}
		else
		{fprintf(stderr,"#Error:Unknown input file format,  %s",filename);exit(0);}
		fclose(fp);
	}
	else
	{fprintf(stderr,"#Error: the file %s could not be opened for reading\n",filename);exit(0);
	}
}//end of init_from_file

int GInitConfig::Output_VTK_header(ofstream *p_fout,int nxx,int nyy,int nzz)
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


