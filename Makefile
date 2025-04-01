CC = gcc
CFLAGS = -march=native -O3
DEBUG_FLAGS = -DDEBUG -g

EXE1 = bin/csyncmer

EXE1_SRCS = src/test.c src/csyncmers.c src/hashing.c src/utils.c

EXE1_LIBS = -lz

EXE1_OBJS = $(EXE1_SRCS:.c=.o)

all: $(EXE1)

debug: CFLAGS += $(DEBUG_FLAGS)
debug: clean $(EXE1)

$(EXE1): $(EXE1_OBJS)
	$(CC) $(CFLAGS) -o $(EXE1) $(EXE1_OBJS) $(EXE1_LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXE1) $(EXE1_OBJS) 

.PHONY: all debug clean