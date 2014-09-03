#!/bin/bash
# This script is intepreted by the Bourne Shell, sh
#
# Documentation for SGE is found in:
# http://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html
#
# Tell SGE which shell to run the job script in rather than depending
# on SGE to try and figure it out.
#$ -S /bin/bash
#
# Export all my environment variables to the job
#$ -V
# Tun the job in the same directory from which you submitted it
#$ -cwd
#
# Give a name to the job
#$ -N APF-STRONG-16
#
# *** Specify the number of cores here 
#$ -pe orte 16

# Specify a time limit for the job
#$ -l h_rt=00:02:00
#
# Join stdout and stderr so they are reported in job output file
#$ -j y
#
# Run on the mpi queue, not more than 2 nodes
#$ -q mpi.q
#
# Specifies the circumstances under which mail is to be sent to the job owner
# defined by -M option. For example, options "bea" cause mail to be sent at the 
# begining, end, and at abort time (if it happens) of the job.
# Option "n" means no mail will be sent.
#$ -m aeb
#
# *** Change to the address you want the notification sent to
#$ -M yourLogin@eng.ucsd.edu
#

# Change to the directory where the job was submitted from
cd $SGE_O_WORKDIR

echo
echo " *** Current working directory"
pwd
echo
echo " *** Compiler"
# Output which  compiler are we using and the environment
mpicc -v
echo
echo " *** Environment"
printenv

echo

echo ">>> Job Starts"
date

# Commands go here

module load rocks-openmpi_openib

mpirun -np 4 ./apf -n 2866 -i 60 -x 1 -y 4

mpirun -np 4 ./apf -n 2866 -i 60 -x 2 -y 2 

mpirun -np 4 ./apf -n 2866 -i 60 -x 4 -y 1


mpirun -np 8 ./apf -n 2866 -i 60 -x 1 -y 8

mpirun -np 8 ./apf -n 2866 -i 60 -x 2 -y 4

mpirun -np 8 ./apf -n 2866 -i 60 -x 4 -y 2

mpirun -np 8 ./apf -n 2866 -i 60 -x 8 -y 1


mpirun -np 16 ./apf -n 2866 -i 60 -x 1 -y 16

mpirun -np 16 ./apf -n 2866 -i 60 -x 2 -y 8

mpirun -np 16 ./apf -n 2866 -i 60 -x 4 -y 4

mpirun -np 16 ./apf -n 2866 -i 60 -x 8 -y 2

mpirun -np 16 ./apf -n 2866 -i 60 -x 16 -y 1

echo "-------------- NO COMM ---------------"
mpirun -np 4 ./apf -n 2866 -i 60 -x 1 -y 4 -k
mpirun -np 4 ./apf -n 2866 -i 60 -x 2 -y 2 -k
mpirun -np 4 ./apf -n 2866 -i 60 -x 4 -y 1 -k
mpirun -np 8 ./apf -n 2866 -i 60 -x 1 -y 8 -k
mpirun -np 8 ./apf -n 2866 -i 60 -x 2 -y 4 -k
mpirun -np 8 ./apf -n 2866 -i 60 -x 4 -y 2 -k
mpirun -np 8 ./apf -n 2866 -i 60 -x 8 -y 1 -k
mpirun -np 16 ./apf -n 2866 -i 60 -x 1 -y 16 -k
mpirun -np 16 ./apf -n 2866 -i 60 -x 2 -y 8 -k
mpirun -np 16 ./apf -n 2866 -i 60 -x 4 -y 4 -k
mpirun -np 16 ./apf -n 2866 -i 60 -x 8 -y 2 -k
mpirun -np 16 ./apf -n 2866 -i 60 -x 16 -y 1 -k

date
echo ">>> Job Ends"
