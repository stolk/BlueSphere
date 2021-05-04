CC=clang-11
#CC=gcc-10

SRC=\
	bluesphere.c \
	pseudorand.c \
	threadpooltask.c \
	threadpool.c \
	

bluesphere: $(SRC) pseudorand.h threadpooltask.h threadpool.h
	$(CC) -Wall -g -O2 -mavx512f -o bluesphere $(SRC) -lm -lpthread


