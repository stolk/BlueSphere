#CC?=clang
CC=gcc

#SANI=-fsanitize=address -fno-omit-frame-pointer

SRC=\
	bluesphere.c \
	pseudorand.c \
	threadpooltask.c \
	threadpool.c \

CFLAGS=\
	-Wall -Wextra -g -O3 -mavx512f $(SANI)

bluesphere: Makefile $(SRC) pseudorand.h threadpooltask.h threadpool.h
	$(CC) $(CFLAGS) -o bluesphere $(SRC) -lm -lpthread

f32tof16: f32tof16.c
	$(CC) $(CFLAGS) -o f32tof16 f32tof16.c -lm

