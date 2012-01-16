CC = gfortran
CFLAGS =
SRC = jacobi1.out gauss-seidel1.out jacobi2.out gauss-seidel2.out cg.out

.SUFFIXES: .f95 .out
.f95.out:
	$(CC) $(CFLAGS) $< -o $@

all: $(SRC)

jacobi2:
	$(CC) $(CFLAGS) -o jacobi2.out jacobi2.f95

gauss-seidel2:
	$(CC) $(CFLAGS) -o gauss-seidel2.out gauss-seidel2.f95

clean:
	rm -f *.out
