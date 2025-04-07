CC = gcc
CFLAGS = -march=native -O3
DEBUG_FLAGS = -DDEBUG -g

EXE1 = bin/csyncmer
EXE2 = bin/syncmer

EXE1_SRCS = src/test.c src/csyncmers.c src/hashing.c src/utils.c
EXE2_SRCS = src/sync.c src/csyncmers.c src/hashing.c src/utils.c src/fasta_reader.c

EXE1_LIBS = -lz
EXE2_LIBS = -lz

EXE1_OBJS = $(EXE1_SRCS:.c=.o)
EXE2_OBJS = $(EXE2_SRCS:.c=.o)

sync: $(EXE2)

all: $(EXE1)

debug: CFLAGS += $(DEBUG_FLAGS)
debug: clean $(EXE1)

$(EXE1): $(EXE1_OBJS)
	$(CC) $(CFLAGS) -o $(EXE1) $(EXE1_OBJS) $(EXE1_LIBS)

$(EXE2): $(EXE2_OBJS)
	$(CC) $(CFLAGS) -o $(EXE2) $(EXE2_OBJS) $(EXE2_LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXE1) $(EXE1_OBJS) $(EXE2) $(EXE2_OBJS)

.PHONY: all debug clean sync