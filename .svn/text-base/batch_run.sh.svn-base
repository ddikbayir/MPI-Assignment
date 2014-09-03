cd /home/rbull/A4/Turnin
module load rocks-openmpi_openib
echo
echo " *** Current working directory"
pwd
echo
echo " *** Compiler"
# Output which  compiler are we using and the environment
gcc -v
echo
echo " *** Environment"
printenv

echo

echo ">>> Job Starts"
date

# Commands go here
./apf -n 200 -i 200
./apf -n 200 -i 200
./apf -n 200 -i 200

./apf -n 400 -i 200
./apf -n 400 -i 200
./apf -n 400 -i 200

./apf -n 500 -i 200
./apf -n 500 -i 200
./apf -n 500 -i 200

mpirun -np 2 ./apf -n 200 -i 200 -x 2
mpirun -np 2 ./apf -n 200 -i 200 -x 2
mpirun -np 2 ./apf -n 200 -i 200 -x 2

mpirun -np 2 ./apf -n 400 -i 200 -x 2
mpirun -np 2 ./apf -n 400 -i 200 -x 2
mpirun -np 2 ./apf -n 400 -i 200 -x 2

mpirun -np 2 ./apf -n 500 -i 200 -x 2
mpirun -np 2 ./apf -n 500 -i 200 -x 2
mpirun -np 2 ./apf -n 500 -i 200 -x 2

mpirun -np 4 ./apf -n 200 -i 200 -x 2 -y 2
mpirun -np 4 ./apf -n 200 -i 200 -x 2 -y 2
mpirun -np 4 ./apf -n 200 -i 200 -x 2 -y 2

mpirun -np 4 ./apf -n 400 -i 200 -x 2 -y 2
mpirun -np 4 ./apf -n 400 -i 200 -x 2 -y 2
mpirun -np 4 ./apf -n 400 -i 200 -x 2 -y 2

mpirun -np 4 ./apf -n 500 -i 200 -x 2 -y 2
mpirun -np 4 ./apf -n 500 -i 200 -x 2 -y 2
mpirun -np 4 ./apf -n 500 -i 200 -x 2 -y 2

mpirun -np 4 ./apf -n 200 -i 200 -y 4
mpirun -np 4 ./apf -n 200 -i 200 -y 4
mpirun -np 4 ./apf -n 200 -i 200 -y 4

mpirun -np 4 ./apf -n 400 -i 200 -y 4
mpirun -np 4 ./apf -n 400 -i 200 -y 4
mpirun -np 4 ./apf -n 400 -i 200 -y 4

mpirun -np 4 ./apf -n 500 -i 200 -y 4
mpirun -np 4 ./apf -n 500 -i 200 -y 4
mpirun -np 4 ./apf -n 500 -i 200 -y 4

mpirun -np 8 ./apf -n 200 -i 200 -y 8
mpirun -np 8 ./apf -n 200 -i 200 -y 8
mpirun -np 8 ./apf -n 200 -i 200 -y 8

mpirun -np 8 ./apf -n 400 -i 200 -y 8
mpirun -np 8 ./apf -n 400 -i 200 -y 8
mpirun -np 8 ./apf -n 400 -i 200 -y 8

mpirun -np 8 ./apf -n 500 -i 200 -y 8
mpirun -np 8 ./apf -n 500 -i 200 -y 8
mpirun -np 8 ./apf -n 500 -i 200 -y 8

mpirun -np 8 ./apf -n 200 -i 200 -x 2 -y 4
mpirun -np 8 ./apf -n 200 -i 200 -x 2 -y 4
mpirun -np 8 ./apf -n 200 -i 200 -x 2 -y 4

mpirun -np 8 ./apf -n 400 -i 200 -x 2 -y 4
mpirun -np 8 ./apf -n 400 -i 200 -x 2 -y 4
mpirun -np 8 ./apf -n 400 -i 200 -x 2 -y 4

mpirun -np 8 ./apf -n 500 -i 200 -x 2 -y 4
mpirun -np 8 ./apf -n 500 -i 200 -x 2 -y 4
mpirun -np 8 ./apf -n 500 -i 200 -x 2 -y 4

mpirun -np 2 ./apf -n 500 -i 200 -x 2 -k
mpirun -np 2 ./apf -n 500 -i 200 -x 2 -k
mpirun -np 2 ./apf -n 500 -i 200 -x 2 -k

mpirun -np 4 ./apf -n 500 -i 200 -x 2 -y 2 -k
mpirun -np 4 ./apf -n 500 -i 200 -x 2 -y 2 -k
mpirun -np 4 ./apf -n 500 -i 200 -x 2 -y 2 -k

mpirun -np 4 ./apf -n 500 -i 200 -y 4 -k
mpirun -np 4 ./apf -n 500 -i 200 -y 4 -k
mpirun -np 4 ./apf -n 500 -i 200 -y 4 -k

mpirun -np 8 ./apf -n 500 -i 200 -y 8 -k
mpirun -np 8 ./apf -n 500 -i 200 -y 8 -k
mpirun -np 8 ./apf -n 500 -i 200 -y 8 -k

mpirun -np 8 ./apf -n 500 -i 200 -x 2 -y 4 -k
mpirun -np 8 ./apf -n 500 -i 200 -x 2 -y 4 -k
mpirun -np 8 ./apf -n 500 -i 200 -x 2 -y 4 -k

