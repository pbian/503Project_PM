SRCS = main_serial.c main_parallel.c
OBJS = $(SRCS:.c=.o)

## turn on when benchmarking, otherwise use -g for debugging
## OPT_FLAG = -O2
OPT_FLAG = -g

CC = gcc

.SUFFIXES: .c .o
.c.o:
	$(CC) -c -o $@ $<
all: particles_serial particles_parallel

particles_parallel: main_parallel.o
	$(CC) -o $@ $^ -lm
particles_serial: main_serial.o
	$(CC) -o $@ $^ -lm
clean:
	/bin/rm -f *.o particles_serial particles_parallel
	/bin/rm -rf particles_serial_program particles_parallel_program
setup:
	mkdir particles_serial_program
	mkdir particles_parallel_program
	mv particles_serial particles_serial_program
	mv particles_parallel particles_parallel_program


