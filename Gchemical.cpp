////////////////////////
// Gchemical.cpp	////
// Y.Gao 05-18-2009	////
////////////////////////

#include "Gchemical.h"
#include <iomanip>
#include "math.h"

Gchemical::Gchemical()
{
	int i;
	ng=nv=nx=ny=nz=0;
	for(i=0;i<10;i++)
		a[i]=0;
	initflag=false;
	option=P234;
}

Gchemical::~Gchemical()
{
}

int Gchemical::Set(int n[],float *aa,chemical_option op)
{
	Set_n(n);
	option=op;
	Set_a(aa);

	initflag=true;

	return 0;
}

int Gchemical::Set_n(int n[])
{
	if(n[0]<1 || n[1]<1 || n[2]<1 || n[3]<1 || n[4]<1)
	{
		cout<<"Invalid dimensions in Gchemical!"<<endl;
		exit(1);
	}
	ng=n[0];
	nv=n[1];
	nx=n[2];
	ny=n[3];
	nz=n[4];

	return 0;
}

int Gchemical::Set_a(float *aa)
{	
	int i;
	if(fabs(aa[0])<1E-5)
	{
		for(i=0;i<10;i++)
			a[i]=aa[i];
	}
	else
	{
		a[0]=aa[0];
			for(i=0;i<10;i++)
				a[i]=aa[i];
	/*	switch(option)
		{
		case P234:
			a[2]=a[0];
			a[3]=-12-3*a[0];
			a[4]=12+2*a[0];
			break;
		case P246:
			a[2]=a[0];
			a[4]=-12-4*a[0];
			a[6]=12+3*a[0];
			break;
		default:
			for(i=0;i<10;i++)
				a[i]=aa[i];
			break;
		}
		*/
	}
	return 0;
}

float Gchemical::Energy(float* eta[],int *gs)
{
	float result=0;
	switch (option)
	{
	case P234:
		result=E_P234(eta,gs);
		break;
	case P246:
		result=E_P246(eta,gs);
		break;
		cout<<"No P246 option for Gchemical!"<<endl;
		exit(1);
	case PREP:
		result=E_PREP(eta,gs);
		break;
	default:
		cout<<"Invalid chemical option!"<<endl;
		exit(1);
	}
	return result;
}

int Gchemical::Potential(float* eta[],float* d_eta[],int *gs,float *conc1,float *conc2)
{
int i,j,k,v1,g1;
long int indx1, indx2;
	switch (option)
	{
	case P234:
		Miu_P234(eta,d_eta,gs);
		break;
	case P246:
		Miu_P246(eta,d_eta,gs,conc1,conc2);
				break;
		cout<<"P246 option for Gchemical!"<<endl;
//		exit(1);
	case PREP:
		Miu_PREP(eta,d_eta,gs);
	//	cout<<"miu of ppt option for Gchemical!"<<endl;
		break;
	default:
		cout<<"Invalid chemical option!"<<endl;
		exit(1);
	}
/*
	char fn[]="miu_v0.vtk";
						for(g1=0;g1<ng;g1++)
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
			fn[5]=(char)((int)'0'+g1*nv+v1);
			ofstream fout(fn,ios::out);
			Output_VTK_header(&fout,nx,ny,nz);
			fout<<fixed;
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);
						fout<<setprecision(6)<<(d_eta[indx1])[indx2]<<endl;
					}
					}
*/
return 0;
}

int Gchemical::UpdateT(float *aa)
{
int i,j,k,v1,g1;
long int indx1, indx2;
		//Miu_P246(eta,d_eta,gs);
	//	a[0]=aa[0];
		a[1]=a[1]+aa[7];
              cout<<"a1="<<a[1]<<endl;
	

return 0;
}

float Gchemical::E_P234(float* eta[],int *gs)
{
	int g1,v1,i,j,k,indx1,indx2;
	float e1,e2,e3,temp;

	float result=0;
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				temp=0;
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])	
						continue;
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						e1=(eta[indx1])[indx2];
						e2=e1*e1;
						e3=e2*e1;
						temp+=e2;
						result+=a[2]/2.0*e2+a[3]/3.0*e3;
					}
					result+=a[4]/4.0*temp*temp;
				}
			}

	return result;
}

float Gchemical::E_P246(float* eta[],int* gs)
{
	float result=0;
	return result;
}

float Gchemical::E_PREP(float* eta[],int* gs)
{
	int i,j,k,g1,v1,g2,v2,indx1,indx1_1,indx2;
	float e1,e2,e3,tp,hs,g;
	float x_Ni,x_Pt,G_B2,G_p,G_total;

	float x1e=0.375,x2e=0.167;
	float ic=0; // |ic|<2

	float dG_1=a[2]+a[4]*2*0.285+a[5]*0.207,
		dG_2=a[3]+a[5]*0.285+a[6]*2*0.207,
		G_00=(a[1]+a[2]*0.285+a[3]*0.207+a[4]*0.285*0.285+a[5]*0.285*0.207+a[6]*0.207*0.207);
	float x1=-(2*dG_1-ic*dG_2)/(2.0*2-ic*ic)/a[0]+0.375,
		x2=-(2*dG_2-ic*dG_1)/(2.0*2-ic*ic)/a[0]+0.167,
		G_p0=G_00-a[0]*((x1-0.375)*(x1-0.375)+(x2-0.167)*(x2-0.167)+ic*(x1-0.375)*(x2-0.167))+0.09*dG_1+(-0.04)*dG_2;

	float result=0;
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);

				x_Ni=eta[index1(0,4)][indx2];
				x_Pt=eta[index1(0,5)][indx2];

				hs=0;
				g=0;
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<4;v1++)
					{
						indx1=index1(g1,v1);
						e1=eta[indx1][indx2];
						e2=e1*e1;
						e3=e2*e1;

						hs+=e3*(e2*6-e1*15+10);
						
						for(g2=0;g2<ng;g2++)
							for(v2=0;v2<4;v2++)
							{
								if(v1==v2)
									continue;
								indx1_1=index1(g2,v2);
								g+=e2*eta[indx1_1][indx2]*eta[indx1_1][indx2]*a[8];
							}
						g+=a[7]*e2*(1-2*e1+e2);
					}
//				G at equil composition
//				G_B2=a[1]+a[2]*0.285+a[3]*0.207+a[4]*0.285*0.285+a[5]*0.285*0.207+a[6]*0.207*0.207;
//				G_p=a[0]*((0.375-x1)*(0.375-x1)+(0.167-x2)*(0.167-x2))+G_p0;

				G_B2=a[1]+a[2]*x_Ni+a[3]*x_Pt+a[4]*x_Ni*x_Ni+a[5]*x_Ni*x_Pt+a[6]*x_Pt*x_Pt;
				G_p=a[0]*((x_Ni-x1)*(x_Ni-x1)+(x_Pt-x2)*(x_Pt-x2))+G_p0;

				G_total=G_B2*(1-hs)+G_p*hs+g;
				result+=G_total;
			}
	
	
	return result;
}

int Gchemical::Miu_P234(float *eta[],float* d_eta[],int *gs)
{
	int g1,v1,i,j,k,indx1,indx2;
	float e1,e2,temp;
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				temp=0;
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])	
						continue;
                                   
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						e1=(eta[indx1])[indx2];
						e2=e1*e1;
						temp+=e2;
					}
                                     
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						e1=(eta[indx1])[indx2];
						(d_eta[v1])[indx2]+=e1*(a[2]+a[3]*e1+a[4]*temp);
					//	(d_eta[v1])[indx2]+=0.05*e1*(2-6*e1+4*temp);
					}
				}
			}
	
	return 0;
}

int Gchemical::Miu_P246(float* eta[],float* d_eta[],int *gs, float *conc1,float *conc2)
{
	int g1,v1,i,j,k,indx1,indx2;
	float e1,e2,temp;
	float x_Ni,x_Pt;
	float T;
	T=a[2];
//        printf("This is miu_246\n");
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++)
        {
				indx2=index2(i,j,k);
          temp=0;
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])	
						continue;
		//			x_Ni=conc1[indx2];
		//			x_Pt=conc2[indx2];
					x_Ni=0.503;
					x_Pt=0.20;
					if(x_Ni<0.5)
					x_Ni=0.5;
               
					for(int v1=0;v1<nv;v1++)
            {
                                                if(v1>9)
                                                       continue;
						indx1=index1(g1,v1);
              e1 = (eta[indx1])[indx2];
              e2=e1*e1;
              temp+=e2;
            }
                     
					for(int v1=0;v1<nv;v1++)
            {
//                                                if(v1>3)
//                                                      continue;
             //     if(v1<4)
		  {
						indx1=index1(g1,v1);
              e1 = (eta[indx1])[indx2];
            /*
              (d_eta[v1])[indx2] += a[6]*e1*e1;
             }
					for(int v1=0;v1<nv;v1++)
            {
						indx1=index1(g1,v1);
              e1 = (eta[indx1])[indx2];
               d_eta[v1][indx2]*=d_eta[v1][indx2];
               d_eta[v1][indx2]*=e1;
               d_eta[v1][indx2]+=e1*(a[2]+a[4]*e1*e1);
               */
            //  (d_eta[v1])[indx2] = e1*(a[2]+a[4]*e1*e1+a[6]*temp*temp);
              //(d_eta[v1])[indx2] = e1*(0.04*(a[2]-120)+a[4]*e1*e1+a[6]*temp*temp);
             // (d_eta[v1])[indx2] = e1*(a[1]+a[4]*e1*e1+a[6]*temp*temp);
             // (d_eta[v1])[indx2] = e1*(a[1]-11.74*e1*e1+17.39*temp*temp);
//              (d_eta[v1])[indx2] = a[0]*e1*(a[1]-13.93*e1*e1+16.30*temp*temp);
//              (d_eta[v1])[indx2] = a[0]*e1*(a[1]*(a[2]-(7541.5-15900*x_Ni+1750*x_Pt))-13.93*e1*e1+16.30*temp*temp);
              (d_eta[v1])[indx2] = a[0]*e1*(a[1]*(a[2]-(7779.5-15900*x_Ni+1750*x_Pt))-1.8213*e1*e1+2.1343*temp*temp);
             // (d_eta[v1])[indx2] = a[0]*e1*(a[1]*(a[2]-(7779.5-15900*eta[10][indx2]+1750*eta[11][indx2]))-1.8213*e1*e1+2.1343*temp*temp);
             // (d_eta[v1])[indx2] = e1*(3.08-12.8*e1*e1+12.6*temp*temp);
            //  (d_eta[v1])[indx2] = e1*(0.1-12.8*e1*e1+12.6*temp*temp);

              // cross-term
              // 
              /*
              if(eta[indx1][indx2]>0.5)
            printf("indx1= %d index2=%d eta=%lf\n",indx1,indx2,eta[indx1][indx2]);
              */
/*
               if(indx1==0 && indx2==68)
{
             printf("eta=%lf\n",eta[0][68]);
             printf("deta=%lf\n",d_eta[0][68]);
             printf("temp=%lf\n",temp);
}
*/
              
							for(int v2=0;v2<nv;v2++)
                {
                  if(v1==v2)
                    continue;
                  (d_eta[v1])[indx2]+=2*e1*pow( (eta[v2])[indx2],2 )*a[8];

                }
             }
	     /*
	     else
	     {
                  (d_eta[v1])[indx2]+=2*(e1-1)*a[7];
							for(int v2=0;v2<nv;v2++)
                {
                  if(v1==v2)
                    continue;
                  (d_eta[v1])[indx2]+=2*e1*pow( (eta[v2])[indx2],2 )*a[8];

                }
            }
	    */
        }//loop of variants
	}
      }//end of space loop
           
	return 0;
}

int Gchemical::Miu_PREP(float* eta[],float* d_eta[],int *gs)
{
	int g1,v1,g2,v2,i,j,k,indx1,indx1_1,indx2;
	float e1,e2,e3,e4,h,dh,hs,dg;
	float h1,dh1,h1s,h2,dh2,h2s;
	float x_Ni,x_Pt,G_B2,G_p;

	float x1e=0.375,x2e=0.167;
	float ic=0; // |ic|<2

	float dG_1=a[2]+a[4]*2*0.285+a[5]*0.207,
		dG_2=a[3]+a[5]*0.285+a[6]*2*0.207,
		G_00=(a[1]+a[2]*0.285+a[3]*0.207+a[4]*0.285*0.285+a[5]*0.285*0.207+a[6]*0.207*0.207);
	float x1=-(2*dG_1-ic*dG_2)/(2.0*2-ic*ic)/a[0]+0.375,
		x2=-(2*dG_2-ic*dG_1)/(2.0*2-ic*ic)/a[0]+0.167,
		G_p0=G_00-a[0]*((x1-0.375)*(x1-0.375)+(x2-0.167)*(x2-0.167)+ic*(x1-0.375)*(x2-0.167))+0.09*dG_1+(-0.04)*dG_2;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
 			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
					continue;

					x_Ni=eta[index1(g1,6)][indx2];
					x_Pt=eta[index1(g1,7)][indx2];
					hs=0;
					h1s=0;
					h2s=0;

					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
						switch(v1)
						{
						case 0:
							
						case 1:
						case 2:
						case 3:
						case 4:
						case 5:
							e1=eta[indx1][indx2];
							e2=e1*e1;
							e3=e2*e1;
							e4=e3*e1;

							h=e3*(10-15*e1+6*e2);
							dh=(e2-2*e3+e4)*30;
							hs+=h;

							h1=e3*(-3*e1+4);
							dh1=12*e2*(-e1+1);
							h1s+=h1;

							h2=e2*(3*e2-8*e1+6);
							dh2=12*e1*(e2-2*e1+1);
							h2s+=h2;

							dg=0;
							for(g2=0;g2<ng;g2++)
								for(v2=0;v2<6;v2++)
								{
									if(v1==v2)
										continue;
									indx1_1=index1(g2,v2);
									dg+=2*e1*eta[indx1_1][indx2]*eta[indx1_1][indx2];
								}
							dg*=a[8];
							dg+=a[7]*(2*e1-6*e2+4*e3);

//							G_B2=a[1]+a[2]*x_Ni+a[3]*x_Pt+a[4]*x_Ni*x_Ni+a[5]*x_Ni*x_Pt+a[6]*x_Pt*x_Pt;
					               // G_B2=a[1]*(2.035*(x_Ni-0.5125)*(x_Ni-0.5125)+0.246*(x_Pt-0.20)*(x_Pt-0.20)-0.0019);
					                G_B2=a[1]*(0.8405*(x_Ni-0.5)*(x_Ni-0.5)+0.123*(x_Pt-0.3063)*(x_Pt-0.3063)-0.0821);
//					                G_B2=a[1]*(2.035*(x_Ni-0.5125)*(x_Ni-0.5125)-0.0019);
//							G_p=a[0]*((x_Ni-x1)*(x_Ni-x1)+(x_Pt-x2)*(x_Pt-x2)+ic*(x_Ni-x1)*(x_Pt-x2))+G_p0;
							G_p=a[0]*(0.8405*(x_Ni-0.53)*(x_Ni-0.53)+0.123*(x_Pt-0.3843)*(x_Pt-0.3843)-0.0845);
	//						G_p=a[0]*(2.035*(x_Ni-0.5425)*(x_Ni-0.5425)-0.0023);
//							G_B2=a[0]*(x_Ni-0.285)*(x_Ni-0.285);
//							G_p=a[0]*(x_Ni-0.375)*(x_Ni-0.375);

//							if(G_B2<0 || G_p<0)
//								cout<<"Need larger G!"<<endl;

							//d_eta[v1][indx2]+=-dh*G_B2+dh*G_p+dg+3.7466E-4*a[3]*e1;
							d_eta[v1][indx2]+=-dh*G_B2+dh*G_p+dg;
							break;
						case 6:	// x_Ni
//							d_eta[v1][indx2]+=(1-hs)*(a[2]+2*a[4]*x_Ni+a[5]*x_Pt)+hs*(2*a[0]*(x_Ni-x1)+a[0]*ic*(x_Pt-x2));
						//	d_eta[v1][indx2]+=(1-hs)*(2*a[1]*2.035*(x_Ni-0.5125))+hs*2.035*(2*a[0]*(x_Ni-0.5425));
							d_eta[v1][indx2]+=(1-hs)*(2*a[1]*0.8405*(x_Ni-0.5))+hs*0.8405*(2*a[0]*(x_Ni-0.53));
							break;
						case 7: // x_Hf
//							d_eta[v1][indx2]+=(1-hs)*(a[3]+a[5]*x_Ni+2*a[6]*x_Pt)+hs*(2*a[0]*(x_Pt-x2)+a[0]*ic*(x_Ni-x1));
							d_eta[v1][indx2]+=(1-hs)*(2*a[1]*0.123*(x_Pt-0.3063))+hs*(2*a[0]*0.123*(x_Pt-0.3843));
							break;
						default:
							cout<<"Invalid variant ID!_Gchemical"<<endl;
							exit(1);
						}
					}
				}
			
			}
	return 0;
}

int Gchemical::Output_VTK_header(ofstream *p_fout,int nxx,int nyy,int nzz)
{
	*p_fout<<"# vtk DataFile Version 2.0"<<endl;
	*p_fout<<"Volume example"<<endl<<"ASCII"<<endl<<"DATASET STRUCTURED_POINTS"<<endl;
	*p_fout<<"DIMENSIONS"<<'\t'<<nxx<<'\t'<<nyy<<'\t'<<nzz<<endl;
	*p_fout<<"DIMENSIONS"<<" "<<nxx<<" "<<nyy<<" "<<nzz<<endl;
	*p_fout<<"ASPECT_RATIO 1 1 1"<<endl;
	*p_fout<<"ORIGIN 0 0 0"<<endl;
	*p_fout<<"POINT_DATA"<<" "<<nxx*nyy*nzz<<endl;
	*p_fout<<"SCALARS volume_scalars float 1"<<endl;
	*p_fout<<"LOOKUP_TABLE default"<<endl;

	return 0;
}

