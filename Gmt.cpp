////////////////////////
// Gmt.cpp	////////////
// Y.Gao 05-20-2009	////
////////////////////////

#include "Gmt.h"
#include <stdio.h>

Gmt::Gmt(char inputfile[])
{	
	Load_Input(inputfile);

	chem.Set(n,chem_a,(chemical_option)chem_op);
	grad.Set(n,kappa,(gradient_option)grad_op);
	elas.Set(n,elas_c,elas_a,rotation,s_app_mag,s_app_inc,e_app_mag,elastic_scale,
		(boundary_condition)bc,(elastic_couple)ec,cd,(fftw_option)fftw_op,
		output_stress,output_interenergy);

	Init_Space();
	Clear_eta(eta_new);
	init.Set(n,(init_option)init_op,eta_old,init_a,input_eta,gs,input_gs,conc1,conc2);
	Set_g();
	Check_gs();
}

Gmt::~Gmt()
{
	Destroy_Space();
}

int Gmt::Load_Input(char inputfile[])
{
	char title[255];
	ifstream fin(inputfile,ios::in);
	if(!fin)
	{
		cerr<<"Input File can not be found!";
		exit(1);
	}
//	Global parameter
	fin>>title;
	fin>>ng>>nv>>nx>>ny>>nz;
	n[0]=ng;
	n[1]=nv;
	n[2]=nx;
	n[3]=ny;
	n[4]=nz;

	cout<<title<<endl;
	cout<<ng<<' '<<nv<<' '<<nx<<' '<<ny<<' '<<nz<<endl;
//	Chemical parameter
	fin>>title;
	fin>>chem_op>>chem_a[0]>>chem_a[1]>>chem_a[2]>>chem_a[3]>>chem_a[4]>>chem_a[5]>>chem_a[6]>>chem_a[7]>>chem_a[8];
	
	cout<<title<<endl;
	cout<<chem_op<<'\t'<<chem_a[0]<<' '<<chem_a[1]<<' '<<chem_a[2]<<' '<<chem_a[3]<<' '<<chem_a[4]<<' '<<chem_a[5]<<' '<<chem_a[6]<<' '<<chem_a[7]<<' '<<chem_a[8]<<endl;
//	Gradient parameter
	fin>>title;
	kappa=new float[nv*nv];
	fin>>grad_op>>kappa[0];
	
	cout<<title<<endl;
	cout<<grad_op<<'\t'<<kappa[0]<<endl;
//	Initial condition
	fin>>title;
	input_gs=new char[255];
	fin>>init_op>>init_a[0]>>init_a[1]>>init_a[2]>>input_eta>>input_gs;

	cout<<title<<endl;
	cout<<init_op<<'\t'<<init_a[0]<<' '<<init_a[1]<<' '<<init_a[2]<<'\t'<<input_eta<<'\t'<<input_gs<<endl;

	if(strcmp(input_gs,"NULL")==0)
	{	
		delete[] input_gs;
		input_gs=NULL;
	}

//	Elastic parameter
	fin>>title;
	fin>>elas_c[0]>>elas_c[1]>>elas_c[2];

	cout<<title<<endl<<elas_c[0]<<' '<<elas_c[1]<<' '<<elas_c[2]<<endl;
	
	fin>>title;
	fin>>elas_a[0]>>elas_a[1]>>elas_a[2]>>elas_a[3];//for monoclinic

//	fin>>elas_a[0]>>elas_a[1]>>elas_a[2];//for orthorhombic
	
	cout<<title<<endl<<elas_a[0]<<' '<<elas_a[1]<<' '<<elas_a[2]<<' '<<elas_a[3]<<endl;//for monoclinic

//	cout<<title<<endl<<elas_a[0]<<' '<<elas_a[1]<<' '<<elas_a[2]<<endl;//for orthorhombic

	fin>>title;
	fin>>s_app_mag>>s_app_inc;

	cout<<title<<endl<<s_app_mag<<'\t'<<s_app_inc<<endl;

	fin>>title;
	fin>>e_app_mag>>e_app_inc;

	cout<<title<<endl<<e_app_mag<<'\t'<<e_app_inc<<endl;

	fin>>title;
	fin>>bc>>ec>>cd;

	cout<<title<<endl<<bc<<'\t'<<ec<<'\t'<<cd<<endl;

	fin>>title;
	fin>>elastic_scale;

	cout<<title<<endl<<elastic_scale<<endl;

	fin>>title;
	fin>>fftw_op;

	cout<<title<<endl<<fftw_op<<endl;
//	Time control
	fin>>title;
	fin>>total_step>>time_step>>output_step;

	cout<<title<<endl<<total_step<<'\t'<<time_step<<'\t'<<output_step<<endl;
//	Random control
	fin>>title;
	fin>>random_coef1>>random_coef2>>random_seed>>random_ave>>random_var>>cut_off_step1>>cut_off_step2>>fluc_space;
	if(nx%fluc_space!=0)
		cout<<"Improper fluc_space!"<<endl;

	cout<<title<<endl<<random_coef1<<' '<<random_coef2<<'\t'<<random_seed<<'\t'<<random_ave<<'\t'<<random_var<<'\t'<<cut_off_step1<<' '<<cut_off_step2<<'\t'<<fluc_space<<endl;
//	Output stress
//        cout<<"Gaussian="<<unirand()<<" "<<unirand()<<endl;
	fin>>title;
	fin>>output_stress;

	cout<<title<<endl<<output_stress<<endl;
//	Output InterEnergy
	fin>>title;
	fin>>output_interenergy;

	cout<<title<<endl<<output_interenergy<<endl;

	
	fin.close();

	return 0;
}

int Gmt::Init_Space()
{
	int g1,v1;
	int mm;
	long int indx1,nxyz=nx*ny*nz;
        char dirname_output[]="output_files"; 
        //initialize concentration
	conc1=new float[nxyz];
	assert(conc1!=NULL);
	conc2=new float[nxyz];
	assert(conc2!=NULL);
	gs=new int[nx*ny*nz];
	assert(gs!=NULL);
	eta_new=new float*[ng*nv];
	assert(eta_new!=NULL);
	eta_old=new float*[ng*nv];
	assert(eta_old!=NULL);
	d_eta=new float*[nv];
	assert(d_eta!=NULL);

	eta_new_k=new fftw_complex*[ng*nv];
	assert(eta_new_k!=NULL);
	eta_old_k=new fftw_complex*[ng*nv];
	assert(eta_old_k!=NULL);
	d_eta_k=new fftw_complex*[ng*nv];
	assert(d_eta_k!=NULL);
	
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			eta_new[indx1]=new float[nxyz];
			assert(eta_new[indx1]!=NULL);
			eta_old[indx1]=new float[nxyz];
			assert(eta_old[indx1]!=NULL);
			if(g1==0)
			{
				d_eta[v1]=new float[nxyz];
				assert(d_eta[v1]!=NULL);
			}

			eta_new_k[indx1]=new fftw_complex[nxyz];
			assert(eta_new_k[indx1]!=NULL);
			eta_old_k[indx1]=new fftw_complex[nxyz];
			assert(eta_old_k[indx1]!=NULL);
			d_eta_k[indx1]=new fftw_complex[nxyz];
			assert(d_eta_k[indx1]!=NULL);
		}
	//	Allocate g,g_sqr
	g_sqr=new float[nxyz];
	assert(g_sqr!=NULL);
	for(mm=0;mm<3;mm++)
	{
		g[mm]=new float[nxyz];
		assert(g[mm]!=NULL);
	}
	outdir=&dirname_output[0];

	return 0;
}

int Gmt::Destroy_Space()
{
	int g1,v1;
	int mm;
	long int indx1;
	
	if(input_gs!=NULL)
		delete[] input_gs;
	delete[] kappa;
	delete[] gs;
	input_gs=NULL;
	kappa=NULL;
	gs=NULL;

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			delete[] eta_new[indx1];
			delete[] eta_old[indx1];
			if(g1==0)
				delete[] d_eta[v1];

			delete[] eta_new_k[indx1];
			delete[] eta_old_k[indx1];
			delete[] d_eta_k[indx1];
		}
	delete[] eta_new;
	delete[] eta_old;
	delete[] d_eta;
	delete[] eta_new_k;
	delete[] eta_old_k;
	delete[] d_eta_k;

	eta_new=NULL;
	eta_old=NULL;
	d_eta=NULL;
	eta_new_k=NULL;
	eta_old_k=NULL;
	d_eta_k=NULL;

//	Deallocate g,g_sqr
	delete[] g_sqr;
	g_sqr=NULL;
	for(mm=0;mm<3;mm++)
	{
		delete[] g[mm];
		g[mm]=NULL;
	}

	return 0;
}

int Gmt::Set_g()
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

int Gmt::Evolution()
{
	int g1,v1,i,j,k;
	long int indx1,indx2,xm,xp,ym,yp,zm,zp;
	float L=5,M1=1,M2=1,M=0;
	float tp;
	
	if(output_stress==1)
	{
//		float s_tp[3][3]={{(float)0.01,0,0},{0,(float)0,0},{0,0,(float)0.01}};		
//		elas.AppliedStrain(s_tp);

		elas.Output_Stress(eta_old,gs);
		cout<<"Output Stress Field!"<<endl;
		elas.Output_TranStrain(eta_old,gs);

		return 1;
	}

	if(output_interenergy==1)
	{
		cout<<"Output Interaction Energy Field!"<<endl;
		elas.Reset_SFTS_MT();
		elas.Output_InterEnergy();

		return 2;
	}

	float **temp,s_tp[3][3]={{(float)0.01,0,0},{0,0,0},{0,0,0}};
	int t=0;
	
//	elas.AppliedStress(s_tp);

	for(t=1;t<=total_step;t++)
	{
		Clear_d_eta();
		
		chem.Potential(eta_old,d_eta,gs,conc1,conc2);
		grad.Potential(eta_old,d_eta,gs);
		elas.Potential(eta_old,d_eta,gs);

		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
			//	tp=0;
				indx2=index2(i,j,k);
		/*	
					xm=index2((i-1+nx)%nx,j,k);
					xp=index2((i+1)%nx,j,k);
					ym=index2(i,(j-1+ny)%ny,k);
					yp=index2(i,(j+1)%ny,k);
					zm=index2(i,j,(k-1+nz)%nz);
					zp=index2(i,j,(k+1)%nz);
                  */        
					for(g1=0;g1<ng;g1++)
					{
						if(g1!=gs[indx2])
							continue;
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
                        //           if (v1==0||v1==2||v1==5||v1>6)
                        //           if (v1==1||v1==2||v1==3||v1>=5)
                              //    if (v1==0||v1==3||v1==4||v1==5||v1==6||v1==7||v1>=10)
                       //    if (v1>3)
                         //             continue;
		/*		    
				     if (eta_new[4][indx2]>0.5)
                                      continue;
			 
				     if (eta_new[5][indx2]>0.5)
                                      continue;
				     if (eta_new[6][indx2]>0.5)
                                      continue;
				     if (eta_new[7][indx2]>0.5)
                                      continue;
				     if (eta_new[8][indx2]>0.5)
                                      continue;
				     if (eta_new[9][indx2]>0.5)
                                      continue;
	*/
			/*  
					
							if(v1>3 && v1<6) // Concentration Field
							{
								if(v1==4)
									M=M1;
								else
									M=M2;

								(eta_new[indx1])[indx2]=(eta_old[indx1])[indx2]
 													+M*((d_eta[v1])[xm]+(d_eta[v1])[xp]+(d_eta[v1])[ym]+(d_eta[v1])[yp]+(d_eta[v1])[zm]+(d_eta[v1])[zp]-6*(d_eta[v1])[indx2])*time_step; // conserve field with constant mobility
								
								if(fabs(random_coef1)>1E-5 && t<=cut_off_step1)
									(eta_new[indx1])[indx2]+=time_step*random_coef1*Gauss_rand();
								
								if(eta_new[indx1][indx2]>10 || eta_new[indx1][indx2]<-10)
								{
									cout<<"Out of range for Concentration!"<<endl;
									cout<<"Eta="<<(eta_new[indx1])[indx2]<<""<<"v="<<indx1<<"indx2="<<indx2<<endl;
									exit(1);
								}
							//	
							//	continue;
							}
							else
                                          //end of concentration               
			*/			
								(eta_new[indx1])[indx2]=(eta_old[indx1])[indx2]
													-L*((d_eta[v1])[indx2])*time_step;
							if(fabs(random_coef1)>1E-5 && t<=cut_off_step1)
							{
								tp=time_step*random_coef1*fabs(Gauss_rand());
							//	tp=time_step*random_coef1*(unirand());
								eta_new[indx1][indx2]+=tp;
								eta_new[index1(g1,4)][indx2]+=tp*0.05;
								eta_new[index1(g1,5)][indx2]+=tp*0.2;
							}
							else if (t<=cut_off_step2)
							{
								tp=time_step*random_coef2*fabs(Gauss_rand());
							//	tp=time_step*random_coef1*(unirand());
								eta_new[indx1][indx2]+=tp;
							
							}
								//eta_new[indx1][indx2]=(eta_new[indx1][indx2]>1.0)?1.0:eta_new[indx1][indx2];
								//eta_new[indx1][indx2]=(eta_new[indx1][indx2]<0.0)?0.0:eta_new[indx1][indx2];
							if(fabs((eta_new[indx1])[indx2])>10000)
							{
								cout<<"Out of range for eta!"<<endl;
								exit(1);
							}
						}
					}
				}
	//	CutLargeEta();

		temp=eta_old;
		eta_old=eta_new;
		eta_new=temp;

		if(t%(total_step/output_step)==0)
		{
			Output_eta(output_step);
			cout<<"t = "<<t<<endl;
		}

	}

	ofstream Eout("Energy.txt",ios::out);
	Eout<<"Chemical free energy: "<<chem.E_PREP(eta_old,gs)<<endl;
	Eout<<"Gradient energy: "<<grad.E_ISO_P7(eta_old,gs)<<endl;

	return 0;
}

int Gmt::Evolution_RK2()
{
	int g1,v1,i,j,k;
	long int indx1,indx2,xm,xp,ym,yp,zm,zp;
	float tp;
	float L=10,M=1;
	int out=1;
	if(output_stress==1)
	{
		float s_tp[3][3]={{(float)1.0,0,0},{0,(float)0,0},{0,0,(float)0.0}};		
 	elas.AppliedStress(s_tp);
//		elas.AppliedStrain(s_tp);

		elas.Output_Stress(eta_old,gs);
		cout<<"Output Stress Field!"<<endl;
		elas.Output_Strain(eta_old,gs,out);

		return 1;
	}

	if(output_interenergy==1)
	{
		cout<<"Output Interaction Energy Field!"<<endl;
		elas.Output_InterEnergy();

		return 2;
	}

//        float **temp,s_tp[3][3]={{(float)1.0,0.0,0},{0.0,0.0,0},{0,0,0}};
	float **temp,e_tp[3][3]={{(float)0.0,0,0},{0,1,0},{0,0,0}};
	int t=0;
	
	elas.AppliedStrain(e_tp);
// 	elas.AppliedStress(s_tp);

	for(t=1;t<=total_step;t++)
	{
		Clear_d_eta();

//	mid value
		chem.Potential(eta_old,d_eta,gs,conc1,conc2);
		grad.Potential(eta_old,d_eta,gs);
		elas.Potential(eta_old,d_eta,gs);

		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);
				/*	
					xm=index2((i-1+nx)%nx,j,k);
					xp=index2((i+1)%nx,j,k);
					ym=index2(i,(j-1+ny)%ny,k);
					yp=index2(i,(j+1)%ny,k);
					zm=index2(i,j,(k-1+nz)%nz);
					zp=index2(i,j,(k+1)%nz);
                             */
					for(g1=0;g1<ng;g1++)
					{
						if(g1!=gs[indx2])
							continue;
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
							
                                if (v1>3)
                                      continue;

		     /* 
				     if (eta_new[4][indx2]>0.5)
                                      continue;
			 
				     if (eta_new[5][indx2]>0.5)
                                      continue;
				     if (eta_new[6][indx2]>0.5)
                                      continue;
				     if (eta_new[7][indx2]>0.5)
                                      continue;
				     if (eta_new[8][indx2]>0.5)
                                      continue;
				     if (eta_new[9][indx2]>0.5)
                                      continue;
                      */
				  //   if (eta_new[10][indx2]>0.515)
                                  //    continue;
				  //   if (eta_new[11][indx2]>0.215)
                                  //    continue;
							/*
							if(v1>3) // Concentration Field
							{
								(eta_new[indx1])[indx2]=(eta_old[indx1])[indx2]
 													+1.0/2*M*((d_eta[v1])[xm]+(d_eta[v1])[xp]+(d_eta[v1])[ym]+(d_eta[v1])[yp]+(d_eta[v1])[zm]+(d_eta[v1])[zp]-6*(d_eta[v1])[indx2])*time_step; // conserve field with constant mobility
								
//								if(fabs(random_coef)>1E-5 && t<=cut_off_step)
//									(eta_new[indx1])[indx2]+=time_step*random_coef*Gauss_rand();
								
								if(eta_new[indx1][indx2]>2 || eta_new[indx1][indx2]<-1)
								{
									cout<<"Out of range for Concentration!"<<endl;
									exit(1);
								}
								continue;
							}
							else // eta Field
							*/
								(eta_new[indx1])[indx2]=(eta_old[indx1])[indx2]
														-1.0/2*L*((d_eta[v1])[indx2])*time_step;
							if(fabs(random_coef1)>1E-5 && t<=cut_off_step1)
							{
								tp=time_step*random_coef1*Gauss_rand();
								eta_new[indx1][indx2]-=tp;
					//			eta_new[index1(g1,4)][indx2]+=tp*0.09;
					//			eta_new[index1(g1,5)][indx2]-=tp*0.04;
							}
							else if(fabs(random_coef2)>1E-5 && t<=cut_off_step2)
							{
								tp=time_step*random_coef2*Gauss_rand();
								eta_new[indx1][indx2]-=tp;
					//			eta_new[index1(g1,4)][indx2]+=tp*0.09;
					//			eta_new[index1(g1,5)][indx2]-=tp*0.04;
							}

							if(fabs((eta_new[indx1])[indx2])>100)
							{
								cout<<"Out of range for eta!"<<endl;
								exit(1);
							}
						}
					}
				}//end of space loop

		temp=eta_old;
		eta_old=eta_new;
		eta_new=temp;

//	final value
/*
		chem.Potential(eta_old,d_eta,gs);
		grad.Potential(eta_old,d_eta,gs);
		elas.Potential(eta_old,d_eta,gs);

		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);
					
					xm=index2((i-1+nx)%nx,j,k);
					xp=index2((i+1)%nx,j,k);
					ym=index2(i,(j-1+ny)%ny,k);
					yp=index2(i,(j+1)%ny,k);
					zm=index2(i,j,(k-1+nz)%nz);
					zp=index2(i,j,(k+1)%nz);

					for(g1=0;g1<ng;g1++)
					{
						if(g1!=gs[indx2])
							continue;
						for(v1=0;v1<nv;v1++)
						{
							indx1=index1(g1,v1);
							*/
							/*
							if(v1>3) // Concentration Field
							{
								(eta_new[indx1])[indx2]=(eta_new[indx1])[indx2]
 													+M*((d_eta[v1])[xm]+(d_eta[v1])[xp]+(d_eta[v1])[ym]+(d_eta[v1])[yp]+(d_eta[v1])[zm]+(d_eta[v1])[zp]-6*(d_eta[v1])[indx2])*time_step; // conserve field with constant mobility
								
//								if(fabs(random_coef)>1E-5 && t<=cut_off_step)
//									(eta_new[indx1])[indx2]+=time_step*random_coef*Gauss_rand();
								
								if(eta_new[indx1][indx2]>2 || eta_new[indx1][indx2]<-1)
								{
									cout<<"Out of range for Concentration!"<<endl;
									exit(1);
								}
								continue;
							}
							else
							*/
							/*
								(eta_new[indx1])[indx2]=(eta_new[indx1])[indx2]
													-L*((d_eta[v1])[indx2])*time_step;
							if(fabs(random_coef1)>1E-5 && t<=cut_off_step1)
							{
								tp=time_step*random_coef1*unirand();
								eta_new[indx1][indx2]+=tp;
							//	eta_new[index1(g1,4)][indx2]+=tp*0.09;//x_Ni
							//	eta_new[index1(g1,5)][indx2]-=tp*0.04;//x_Pt
							}
								
							if(fabs((eta_new[indx1])[indx2])>100)
							{
								cout<<"Out of range for eta!"<<endl;
								exit(1);
							}
						}
					}
				}

		temp=eta_old;
		eta_old=eta_new;
		eta_new=temp;
*/


		CutLargeEta();
		if(t%(total_step/output_step)==0)
		{
	//	       elas.UpdateStrain(e_tp);
      	//	      elas.UpdateStress(s_tp);
//		       chem.UpdateT(chem_a);
		  if (t<=total_step/2)
		      { chem.UpdateT(chem_a);}
      		 //    elas.UpdateStress(s_tp);
		       else
		      { 
		      chem.UpdateT2(chem_a);
      		//      elas.UpdateStress2(s_tp);
		       }

			Output_eta(output_step);
	//	elas.Output_Stress(eta_old,gs);
	//	        elas.Output_Strain(eta_old,gs,output_step);
			cout<<"t = "<<t<<endl;
		}

	}

	ofstream Eout("Energy.txt",ios::out);
//	Eout<<"Chemical free energy: "<<chem.E_PREP(eta_old,gs)<<endl;
//	Eout<<"Gradient energy: "<<grad.E_ISO_P7(eta_old,gs)<<endl;

	return 0;
}

int Gmt::Evolution_Spectrum()
{
	int g1,v1,i,j,k;
	long int indx1,indx2,indx22;
	float L=5,M1=1,M2=1,M=0;

	float kap=0;
	int vrand=0;
	float tp=0,tp1=1.1;
	fftw_complex **temp;

	int ii,jj,kk;
	
	Gfftw3 gfft;
	gfft.Set(n,(fftw_option)fftw_op);

	if(output_stress==1)
	{
//		float s_tp[3][3]={{(float)0.01,0,0},{0,(float)0,0},{0,0,(float)0.01}};
//		elas.AppliedStrain(s_tp);

		elas.Output_Stress(eta_old,gs);
		cout<<"Output Stress Field!"<<endl;
		elas.Output_TranStrain(eta_old,gs);
	//	elas.Output_Strain(eta_old,gs);

		return 1;
	}

	if(output_interenergy==1)
	{
		cout<<"Output Interaction Energy Field!"<<endl;

		elas.Reset_SFTS_MT(); // Interaction Energy for MT in NiTiPt
		elas.Output_InterEnergy();

		return 2;
	}

	if(ng!=1)
	{
		cout<<"Invalid Grain number for Spectrum method!"<<endl;
		exit(1);
	}

	float s_tp[3][3]={{0,0,0},{0,0,0},{0,0,0}};
	
//	if(fabs(s_app_mag)>0.00001)
//		elas.AppliedStress(s_tp);

	int t=0;
	
	for(t=1;t<=total_step;t++)
	{
		Clear_d_eta();
		Clear_d_eta_k();
		
		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			//for(v1=0;v1<12;v1++)//Taiwu modified
			{
                              //   if(v1>2)
                              //    continue;
				indx1=index1(g1,v1);
				gfft.FFTW_3D(eta_old[indx1],eta_old_k[indx1],FORWARD);
			}

		chem.Potential(eta_old,d_eta,gs,conc1,conc2);
		elas.Potential(eta_old,d_eta,gs);

		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			//for(v1=0;v1<12;v1++)
			{
                                // if(v1>2)
                                //  continue;
				indx1=index1(g1,v1);
				gfft.FFTW_3D(d_eta[indx1],d_eta_k[indx1],FORWARD);

				for(i=0;i<nx;i++)
					for(j=0;j<ny;j++)
						for(k=0;k<nz;k++)
						{
							indx2=index2(i,j,k);

							if(v1>5) // concentration
							{
								if(v1==6)
								{
									M=M1;
									kap=10*kappa[0];
								}
								else
								{
									M=M2;
									kap=5*kappa[0];
								}

								eta_new_k[indx1][indx2][0]=(eta_old_k[indx1][indx2][0]*(1-0*M*time_step*g_sqr[indx2]*g_sqr[indx2]*kap)-M*time_step*d_eta_k[indx1][indx2][0]*g_sqr[indx2])
															/(1+1*M*time_step*g_sqr[indx2]*g_sqr[indx2]*kap);
								eta_new_k[indx1][indx2][1]=(eta_old_k[indx1][indx2][1]*(1-0*M*time_step*g_sqr[indx2]*g_sqr[indx2]*kap)-M*time_step*d_eta_k[indx1][indx2][1]*g_sqr[indx2])
															/(1+1*M*time_step*g_sqr[indx2]*g_sqr[indx2]*kap);
							}
							else // eta
							
							{
								eta_new_k[indx1][indx2][0]=(eta_old_k[indx1][indx2][0]*(1-0*L*time_step*g_sqr[indx2]*kappa[0])-L*time_step*d_eta_k[indx1][indx2][0])
															/(1+1*L*time_step*g_sqr[indx2]*kappa[0]);
								eta_new_k[indx1][indx2][1]=(eta_old_k[indx1][indx2][1]*(1-0*L*time_step*g_sqr[indx2]*kappa[0])-L*time_step*d_eta_k[indx1][indx2][1])
															/(1+1*L*time_step*g_sqr[indx2]*kappa[0]);
							}
						}
			}

		temp=eta_old_k;
		eta_old_k=eta_new_k;
		eta_new_k=temp;

		for(g1=0;g1<ng;g1++)
			for(v1=0;v1<nv;v1++)
			{
				indx1=index1(g1,v1);
				gfft.FFTW_3D(eta_old_k[indx1],eta_old[indx1],BACKWARD);
			}

		// Random fluctuation in normal space
		// Random fluctuation in normal space
		LangevinNoise(t);
	//	
/*	
							if(fabs(random_coef1)>1E-5 && t<=cut_off_step1)
							{
								tp=time_step*random_coef1*Gauss_rand();
								eta_old[indx1][indx2]+=tp;
								eta_old[index1(g1,4)][indx2]+=tp*0.04;
								eta_old[index1(g1,5)][indx2]+=tp*0.09;
							}
*/					
		CutLargeEta();

/*				for(i=0;i<nx;i++)
					for(j=0;j<ny;j++)
						for(k=0;k<nz;k++)
						{
                eta_old[index1(0,3)][indx2]=eta_old[index1(0,0)][indx2]+eta_old[index1(0,1)][indx2]+eta_old[index1(0,2)][indx2];
                        }
*/
		if(t%(total_step/output_step)==0)
		{
			cout<<"t = "<<t<<endl;
			Output_eta(output_step);
		}

		if(fabs(eta_old[indx1][index2(nx/2,ny/2,nz/2)])>100)
		{
			cout<<"Out of range for spectrum method!"<<endl;
			exit(1);
		}

	}

	ofstream Eout("Energy.txt",ios::out);
	Eout<<"Chemical free energy: "<<chem.Energy(eta_old,gs)<<endl;
	Eout<<"Gradient energy: "<<grad.Energy(eta_old,gs)<<endl;

	return 0;
}

int Gmt::Output_eta(int out)
{
	int g1,v1,i,j,k;
	long int indx1,indx2,indx3;
	float max=-100,min=100,fraction;
	float xNi_ave=0,xPt_ave=0,xNi_max=0,xNi_min=1,xPt_max=0,xPt_min=1;
	float eta_ave[4]={0};
	float eta_total[nx*ny*nz];
	int vID;
	float vmax,temp;
	
	static int num=0;
	char f0[20],f1[20],f2[20],f3[20];
	char f4[20],f5[20],f6[20],f7[20];
	char f8[20],f9[20],f10[20],f11[20];
        char fppt[20],ft[20];
	float **eta_out=new float *[nv];
	assert(eta_out!=NULL);
	for(int i=0;i<nv;i++)
	{
		eta_out[i]=new float [(nx+2)*(ny+2)*(nz+2)];
		assert(eta_out[i]!=NULL);
	}
	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(k=0;k<nz+2;k++)
				for(j=0;j<ny+2;j++)
					for(i=0;i<nx+2;i++)
					{
						indx3=index3(i,j,k);
                            (eta_out[indx1])[indx3]=0;    
				}
}

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
			indx1=index1(g1,v1);
			for(k=1;k<nz+1;k++)
				for(j=1;j<ny+1;j++)
					for(i=1;i<nx+1;i++)
					{
						indx3=index3(i,j,k);
						indx2=index2(i-1,j-1,k-1);
                                    (eta_out[indx1])[indx3]=(eta_old[indx1])[indx2];
				    }
}

//	sprintf(f0,"ther_detw_v1_%02d.vtk",num);
//	sprintf(f1,"ther_detw_v2_%02d.vtk",num);
//	sprintf(f2,"ther_detw_v3_%02d.vtk",num);
//	sprintf(f3,"ther_detw_v4_%02d.vtk",num);
//	sprintf(f0,"test_detwin_T20_v1_%02d.vtk",num);
//	sprintf(f1,"test_detwin_T20_v2_%02d.vtk",num);
//	sprintf(f0,"tetra_v1_3f_%02d.vtk",num);
//	sprintf(f1,"tetra_v2_3f_%02d.vtk",num);
//	sprintf(f2,"tetra_v3_3f_%02d.vtk",num);

/*
	sprintf(f0,"ther_c1_T7r_sc_v1_%02d.vtk",num);
	sprintf(f1,"ther_c1_T7r_sc_v2_%02d.vtk",num);
	sprintf(f2,"ther_c1_T7r_sc_v3_%02d.vtk",num);
	sprintf(f3,"ther_c1_T7r_sc_v4_%02d.vtk",num);
*/

        sprintf(f0,"MT_v1_%02d.vtk",num);
	sprintf(f1,"MT_v2_%02d.vtk",num);
	sprintf(f2,"MT_v3_%02d.vtk",num);
	sprintf(f3,"MT_v4_%02d.vtk",num);

/*
        sprintf(f0,"T150_e12s_v1_%02d.vtk",num);
	sprintf(f1,"T150_e12s_v2_%02d.vtk",num);
	sprintf(f2,"T150_e12s_v3_%02d.vtk",num);
	sprintf(f3,"T150_e12s_v4_%02d.vtk",num);
*/
/*

	sprintf(f4,"T150_e12s_v5_%02d.vtk",num);
	sprintf(f5,"T150_e12s_v6_%02d.vtk",num);
	sprintf(f6,"T150_e12s_v7_%02d.vtk",num);
	sprintf(f7,"T150_e12s_v8_%02d.vtk",num);
	sprintf(f8,"T150_e12s_v9_%02d.vtk",num);
	sprintf(f9,"T150_e12s_v10_%02d.vtk",num);
	sprintf(f10,"T150_e12s_v11_%02d.vtk",num);
	sprintf(f11,"T150_e12s_v12_%02d.vtk",num);
*/
//
/*
	sprintf(f0,"input_v1_part5_%02d.vtk",num);
	sprintf(f1,"input_v2_part5_%02d.vtk",num);
	sprintf(f2,"input_v3_part5_%02d.vtk",num);
	sprintf(f3,"input_v4_part5_%02d.vtk",num);

	sprintf(f4,"input_v5_part5_%02d.vtk",num);
	sprintf(f5,"input_v6_part5_%02d.vtk",num);
//
        sprintf(f6,"input_Ni_part5_%02d.vtk",num);
        sprintf(f7,"input_Hf_part5_%02d.vtk",num);
*/
//	sprintf(f7,"Hphase_two3_Hf_%02d.vtk",num);
	
/*	sprintf(f8,"M_pure_v9_%02d.vtk",num);
	sprintf(f9,"M_pure_v10_%02d.vtk",num);
	sprintf(f10,"M_pure_v11_%02d.vtk",num);
	sprintf(f11,"M_pure_v12_%02d.vtk",num);
	*/
//	sprintf(ft,"test_detwin_T20_total_%02d.vtk",num);
//	sprintf(ft,"tetra_total3f_%02d.vtk",num);

//	sprintf(filename,"Hphase_p3_v1_02.vtk");
//	sprintf(ft,"ther_l2_T2_total_%02d.vtk",num);
	//sprintf(ft,"Mart_parnew5_total_%02d.vtk",num);
  //      ft[7]=char((int)ft[7]+num);
//	char fn[20]="eta00.vtk";
//	fn[4]=f0[4];
        num++;
//	ofstream fout(fn,ios::out);
	ofstream v0_out(f0,ios::out);
	ofstream v1_out(f1,ios::out);
	ofstream v2_out(f2,ios::out);
	ofstream v3_out(f3,ios::out);
/*	ofstream v4_out(f4,ios::out);
	ofstream v5_out(f5,ios::out);
	ofstream v6_out(f6,ios::out);
	ofstream v7_out(f7,ios::out);
*/
/*
	ofstream v8_out(f8,ios::out);
	ofstream v9_out(f9,ios::out);
	ofstream v10_out(f10,ios::out);
	ofstream v11_out(f11,ios::out);
*/
/* 
	ofstream v4_out(f4,ios::out);
	ofstream v5_out(f5,ios::out);

        ofstream v6_out(f6,ios::out);
	ofstream v7_out(f7,ios::out);
*/
/*
	ofstream v8_out(f8,ios::out);
	ofstream v9_out(f9,ios::out);
	ofstream v10_out(f10,ios::out);
	ofstream v11_out(f11,ios::out);
*/
//	ofstream ppt_out(fppt,ios::out);
//	ofstream ft_out(ft,ios::out);

//	ofstream c_out("Concentrations.txt",ios::app);

//	Output_VTK_header(&fout,nx,ny,nz*nv);
	Output_VTK_header(&v0_out,nx+2,ny+2,nz+2);
	Output_VTK_header(&v1_out,nx+2,ny+2,nz+2);
	Output_VTK_header(&v2_out,nx+2,ny+2,nz+2);
	Output_VTK_header(&v3_out,nx+2,ny+2,nz+2);
/*	Output_VTK_header(&v4_out,nx,ny,nz);	
	Output_VTK_header(&v5_out,nx,ny,nz);

        Output_VTK_header(&v6_out,nx,ny,nz);
	Output_VTK_header(&v7_out,nx,ny,nz);
*/
/*
	Output_VTK_header(&v8_out,nx,ny,nz);
	Output_VTK_header(&v9_out,nx,ny,nz);
	Output_VTK_header(&v10_out,nx,ny,nz);
	Output_VTK_header(&v11_out,nx,ny,nz);
*/
//	Output_VTK_header(&ft_out,nx,ny,nz);

	for(g1=0;g1<ng;g1++)
		for(v1=0;v1<nv;v1++)
		{
		if(v1>11)
		continue;
			indx1=index1(g1,v1);
			for(k=0;k<nz;k++)
				for(j=0;j<ny;j++)
					for(i=0;i<nx;i++)
					{
						indx2=index2(i,j,k);

				//	fout<<(eta_old[indx1])[indx2]<<endl;
						
						max=(eta_old[indx1])[indx2]>max?(eta_old[indx1])[indx2]:max;
						min=(eta_old[indx1])[indx2]<min?(eta_old[indx1])[indx2]:min;
					}
		}

	for(k=0;k<nz+2;k++)
		for(j=0;j<ny+2;j++)
			for(i=0;i<nx+2;i++)
			{
                        
				indx3=index3(i,j,k);
		//		vID=-1;
		//		vmax=0;

				for(g1=0;g1<ng;g1++)
				{
					if(g1!=gs[indx2])
						continue;
						/*
					for(v1=0;v1<nv;v1++)
					{
					
						indx1=index1(g1,v1);
						temp=(eta_old[indx1])[indx2];
						if(temp>vmax)
						{
							vmax=temp;
							vID=v1;
						}
					}
					*/
				}
//                      eta_total[indx2]=eta_old[index1(0,0)][indx2]+eta_old[index1(0,1)][indx2]*2+ \
//                                          eta_old[index1(0,2)][indx2]*3+eta_old[index1(0,3)][indx2]*4;
//                      eta_total[indx2]=eta_old[index1(0,0)][indx2]+eta_old[index1(0,1)][indx2]*2+eta_old[index1(0,2)][indx2]*3;
//                                          eta_old[index1(0,2)][indx2]*3+eta_old[index1(0,3)][indx2]*4;
/*
					   eta_old[index1(0,4)][indx2]*5+eta_old[index1(0,5)][indx2]*6+ \
                                           eta_old[index1(0,6)][indx2]*7+eta_old[index1(0,7)][indx2]*8+ \
					   eta_old[index1(0,8)][indx2]*9+eta_old[index1(0,9)][indx2]*10+ \
					   eta_old[index1(0,10)][indx2]*11+eta_old[index1(0,11)][indx2]*12;
*/
//		if(fabs(vmax)<0.5)
		//			vID=-1;
					
			/*		
				v0_out<<eta_old[index1(0,0)][indx2]<<endl;
				v1_out<<eta_old[index1(0,1)][indx2]<<endl;
				v2_out<<eta_old[index1(0,2)][indx2]<<endl;
				v3_out<<eta_old[index1(0,3)][indx2]<<endl;
			
				v4_out<<eta_old[index1(0,4)][indx2]<<endl;//x_Ni
				v5_out<<eta_old[index1(0,5)][indx2]<<endl;//x_Hf
			*/
			/*	v6_out<<eta_old[index1(0,6)][indx2]<<endl;
				v7_out<<eta_old[index1(0,7)][indx2]<<endl;
			*/
					
				v0_out<<eta_out[index1(0,0)][indx3]<<endl;
				v1_out<<eta_out[index1(0,1)][indx3]<<endl;
				v2_out<<eta_out[index1(0,2)][indx3]<<endl;
				v3_out<<eta_out[index1(0,3)][indx3]<<endl;
   	/*	
				v4_out<<eta_old[index1(0,4)][indx2]<<endl;
				v5_out<<eta_old[index1(0,5)][indx2]<<endl;

                                v6_out<<eta_old[index1(0,6)][indx2]<<endl;//x_Ni
//                                v6_out<<conc1[indx2]<<endl;//x_Ni
				v7_out<<eta_old[index1(0,7)][indx2]<<endl;//x_Hf
*/
/*				
		
				v8_out<<eta_old[index1(0,8)][indx2]<<endl;
				v9_out<<eta_old[index1(0,9)][indx2]<<endl;
				v10_out<<eta_old[index1(0,10)][indx2]<<endl;
				v11_out<<eta_old[index1(0,11)][indx2]<<endl;
   */                        
//                          ft_out<<eta_total[indx2]<<endl;

			}//end of space loop`


//	ft_out.flush();
//	ft_out.close();

	
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
						fraction+=(eta_old[indx1])[indx2];
					}
			cout<<fraction/nx/ny/nz<<endl;
		}	
	

	v0_out.flush();
	v0_out.close();
	v1_out.flush();
	v1_out.close();
	v2_out.flush();
	v2_out.close();
	v3_out.flush();
	v3_out.close();
/*	v4_out.flush();
	v4_out.close();

	v5_out.flush();
	v5_out.close();
	
	v6_out.flush();
	v6_out.close();
	v7_out.flush();
	v7_out.close();	
*/
/*
	v8_out.flush();
	v8_out.close();
	v9_out.flush();
	v9_out.close();
	v10_out.flush();
	v10_out.close();
	v11_out.flush();
	v11_out.close();
*/
//        ft_out.flush();
//	ft_out.close();

//	xNi_ave/=nx*ny*nz;
//	xPt_ave/=nx*ny*nz;

	cout<<"max="<<max<<'\t'<<"min="<<min<<endl;

	return 0;
}

int Gmt::Output_VTK_header(ofstream *p_fout,int nxx,int nyy,int nzz)
{
	*p_fout<<"# vtk DataFile Version 2.0"<<endl;
	*p_fout<<"Volume example"<<endl<<"ASCII"<<endl<<"DATASET STRUCTURED_POINTS"<<endl;
	*p_fout<<"DIMENSIONS"<<" "<<nxx<<" "<<nyy<<" "<<nzz<<endl;
	*p_fout<<"ASPECT_RATIO 1 1 1"<<endl;
	*p_fout<<"ORIGIN 0 0 0"<<endl;
	*p_fout<<"POINT_DATA"<<" "<<(nxx)*(nyy)*(nzz)<<endl;
	*p_fout<<"SCALARS volume_scalars float 1"<<endl;
	*p_fout<<"LOOKUP_TABLE default"<<endl;

	return 0;
}

int Gmt::Output_gs()
{
	int i,j,k;
	long indx2;
	ofstream fout("gs.dat",ios::out);
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				fout<<gs[indx2]<<'\t';
			}
	fout.flush();
	fout.close();
	return 0;
}

int Gmt::Clear_d_eta()
{
	int v1,i,j,k;
	long int indx2;
	for(v1=0;v1<nv;v1++)
		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
				for(k=0;k<nz;k++)
				{
					indx2=index2(i,j,k);
					d_eta[v1][indx2]=0;
				}
	return 0;
}

int Gmt::Clear_d_eta_k()
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
						d_eta_k[indx1][indx2][0]=0;
						d_eta_k[indx1][indx2][1]=0;
					}
		}
	return 0;
}

int Gmt::Clear_eta(float *eta[])
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
						(eta[indx1])[indx2]=0;
					}
		}
	return 0;
}

bool Gmt::Check_gs()
{
	int i,j,k;
	long int indx2;
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				if(gs[indx2]>=ng)
				{
					cout<<"Invalid grain shape distribution!"<<endl;
					exit(1);
				}
			}
	return true;
}

int Gmt::LangevinNoise(int time)
{
	int i,j,k,v1;
	int ii,jj,kk;
	long int indx1,indx2,indx22;

	float rc=0;

	float tp=0;

	if(time<=cut_off_step1)
		rc=random_coef1;
	else if(time<=cut_off_step2)
		rc=random_coef2;
	else
		return 1;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				if(fabs(rc)>1E-5 && time<=cut_off_step2)
				{
			//		if(i%fluc_space==0 && j%fluc_space==0 && k%fluc_space==0)
			//		{
						for(v1=0;v1<nv;v1++)
						{
                                  //            if(v1==4||v1==5)
				//	      {
				//	      rc=0.04*rc;
				//	      }
                                              //         continue;
							tp=time_step*rc*Gauss_rand();
					//		tp=time_step*rc*(unirand()-0.5);
					 //    if(v1==0)
					 //    cout<<"random noise="<<tp<<endl;
							
							if(nz>2) // 3D
							{
							
							//	for(ii=0;ii<fluc_space;ii++)
							//		for(jj=0;jj<fluc_space;jj++)
							//			for(kk=0;kk<fluc_space;kk++)
							//			{
							//				indx22=index2(i+ii,j+jj,k+kk);
											indx1=index1(0,v1);
										if(v1==6||v1==7)
									//	if(v1==6)
										{
											eta_old[index1(0,6)][indx2]+=0.04*tp;
											eta_old[index1(0,7)][indx2]+=0.04*tp;
											}
											else{
											eta_old[indx1][indx2]+=tp;
											}
							//			}

							}
							else // 2D: nz=1
							{
								for(ii=0;ii<fluc_space;ii++)
									for(jj=0;jj<fluc_space;jj++)
									{
										indx22=index2(i+ii,j+jj,0);
										indx1=index1(0,v1);
										eta_old[indx1][indx22]+=tp;
										eta_old[indx1][index2(i+ii,j+jj,1)]=eta_old[indx1][index2(i+ii,j+jj,0)];
									}
							}
						}
				//	}
				}
			}

	return 0;
}


int Gmt::CutLargeEta()
{
	int i,j,k,g1,v1;
	long int indx1,indx2;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			for(k=0;k<nz;k++)
			{
				indx2=index2(i,j,k);
				for(g1=0;g1<ng;g1++)
					for(v1=0;v1<nv;v1++)
					{
						indx1=index1(g1,v1);
         /*                          
						if(v1==6)
						{
						if(eta_old[indx1][indx2]<0.498)
							eta_old[indx1][indx2]=0.498;
						
						}
					else if(v1==7)
						{
						if(eta_old[indx1][indx2]<0.19)
							eta_old[indx1][indx2]=0.19;
						
						}
						
					//		continue;
                                                else
			*/			{
						if(eta_old[indx1][indx2]<0)
							eta_old[indx1][indx2]=0;
						if(eta_old[indx1][indx2]>1)
							eta_old[indx1][indx2]=1;
					        }
					}
			}

	return 0;
}
 double Gmt::unirand()
{
srand(0);
return ((double)rand())/RAND_MAX;
}
