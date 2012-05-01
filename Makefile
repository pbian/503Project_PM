SRCS = main_serial.c main_parallel_MPI.c main_serial_2D.c main_parallel_MPI_2D.c
OBJS = $(SRCS:.c=.o)

## turn on when benchmarking, otherwise use -g for debugging
## OPT_FLAG = -O2
OPT_FLAG = -g -mpitrace

CC = openmpicc

.SUFFIXES: .c .o
.c.o:
	$(CC) -c -o $@ $<
all_1D: particles_serial particles_parallel
all_2D: particles_serial_2D particles_parallel_2D

particles_parallel: main_parallel_MPI.o
	$(CC) -o $@ $^ -lm
particles_serial: main_serial.o
	$(CC) -o $@ $^ -lm
particles_serial_2D: main_serial_2D.o
	$(CC) -o $@ $^ -lm
particles_parallel_2D: main_parallel_MPI_2D.o
	$(CC) -o $@ $^ -lm
clean:
	/bin/rm -f *.o particles_serial particles_parallel particles_serial_2D particles_parallel_2D
	/bin/rm -rf particles_serial_program particles_parallel_program particles_serial_2D_program particles_parallel_2D_program
setup_2D:
	mkdir particles_parallel_2D_program
	mkdir particles_serial_2D_program
	mv particles_parallel_2D particles_parallel_2D_program
	mv particles_serial_2D particles_serial_2D_program
setup_1D:
	mkdir particles_serial_program
	mkdir particles_parallel_program
	mv particles_serial particles_serial_program
	mv particles_parallel particles_parallel_program

