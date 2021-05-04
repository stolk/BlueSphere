CC=clang-11
#CC=gcc-10

#SANI=-fsanitize=address -fno-omit-frame-pointer

SRC=\
	bluesphere.c \
	pseudorand.c \
	threadpooltask.c \
	threadpool.c \

CFLAGS=\
	-Wall -g -O2 -mavx512f $(SANI)

bluesphere: Makefile $(SRC) pseudorand.h threadpooltask.h threadpool.h
	$(CC) $(CFLAGS) -o bluesphere $(SRC) -lm -lpthread

f32tof16: f32tof16.c
	$(CC) $(CFLAGS) -o f32tof16 f32tof16.c -lm

