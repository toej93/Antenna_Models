#/bin/bash
#PBS -N Antenna_models
#PBS -l nodes=1:ppn=2
#PBS -j oe
#PBS -A PAS0654
#PBS -m e
#PBS -l mem=64000MB
#PBS -l walltime=01:00:00

source /users/PCON0003/cond0068/.bash_profile_owens


cd $PBS_O_WORKDIR
cd /users/PCON0003/cond0068/ARA/Antenna_Models

root -l ariannaHeff.C &
#root -l Converter_WIPLD_to_ARA.C &

wait
