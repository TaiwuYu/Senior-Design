[number_of_grain,variants,nx,ny,nz]
1 4 128 128 128

[chemical_option,a0,a1,a2,a3,a4,a5,a6,a7,a8]
1 32 0.0094 125  0 -11.74 0 17.39 0.0 0.0

[gradient_option,kappa]
0  6

[initial_condition(0=HOMO,1=SINGLE,2=TWIN_13,3=USER,4=LOAD,5=GRAIN_3,6=GRAIN_4,7=DOUBLE,8=TWIN_12,9=TWIN_159,10=LAM_123,11=LOAD_3264C,12=TWO,13=HPHASE,14=USER_SET,15=RANDOM_PARTICLE),init_a[3],input_file]
3 5 40 5  eta_input.vtk NULL

[elastic_constant:c11,c12,c44:1.83,1.46,0.46]
2.00 1.36 0.49

[lattice_mismatch(4):B2_R,alfa,beta,gamma:0.99978,0.0.99978,0.00603,0.00603;B2_B19_prime,gamma,alfa,epsilon,delta:0.9892,1.0198,-0.0654,0.0908]
0.9892 1.0198 -0.0654 0.0908

[applied_stress_magnitude,applied_stress_increment]
100.0 0.0

[applied_strain_magnitude,applied_strain_increment]
0.0 0.0

[boundary_condition(FIX,RELAX),elastic_couple(PHI1,PHI2,PHI23,PHI345),config-dept]
0	1	3

[elastic_energy_scale:35.2]
1000

[fftw3_option:ESTIMATE,MEASURE,PATIENT,EXHAUSTIVE]
1

[total_step,time_step,output_step]
4000	0.0025  40

[coefficient(1/2),random_seed,average,variation,cut_off_step(1/2),fluc_space]
3.0  0.00 130   0  5 500 4000 1 

[Output_stress/strain]
0

[Output_InterEnergy]
0

