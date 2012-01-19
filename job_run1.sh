# nqs_mpi_fujitsu.sh
# @$-eo
#
# @$-q gh10034
# @$-g gh10034
# @$-lP 1
# @$-lp 16
# @$-lm 1800mb
# @$-cp 0:10:00

set -x

cd $QSUB_WORKDIR
./jacobi1.out
./jacobi1.out
./jacobi1.out
./gauss-seidel1.out
./gauss-seidel1.out
./gauss-seidel1.out
