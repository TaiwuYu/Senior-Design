#PBS -A PAS0971	
#PBS -N Hphase_rand2
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -S /bin/bash
#PBS -m abe

#Job script for owens
module load intel
#module load gnu/4.8.5
module load fftw3/3.3.8

cd $PBS_O_WORKDIR
#export OMP_NUM_THREADS=28
./Hphase_model.exe &> Hphase_rand2.out
