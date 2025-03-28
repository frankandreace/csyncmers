# all:
# 	gcc -march=native -Ofast -o syncmer syncmer.c

CC = gcc
CFLAGS1 = -march=native -O3
CFLAGS2 = -march=native -O3

EXE1 = bin/syncmer_tree
EXE2 = bin/csyncmer

EXE1_SRCS = src/syncmer.c
EXE2_SRCS = src/csyncmer.c src/closed_syncmers_syng.c src/utils.c src/fasta_reader.c

EXE2_LIBS = -lz

EXE1_OBJS = $(EXE1_SRCS:.c=.o)
EXE2_OBJS = $(EXE2_SRCS:.c=.o)

all: $(EXE1) $(EXE2)

$(EXE1): $(EXE1_OBJS)
	$(CC) $(CFLAGS1) -o $(EXE1) $(EXE1_OBJS) 

$(EXE2): $(EXE2_OBJS)
	$(CC) $(CFLAGS2) -o $(EXE2) $(EXE2_OBJS) $(EXE2_LIBS)

%.o: %.c
	$(CC) $(CFLAGS1) -c $< -o $@

clean:
	rm -f $(EXE1) $(EXE2) $(EXE1_OBJS) $(EXE2_OBJS)
