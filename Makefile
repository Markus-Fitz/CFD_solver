# Minimal Makefile for a simple 2D CFD starter (C, no GPU)
CC = gcc
CFLAGS = -std=c11 -O2 -Iinclude -Wall -Wextra
LDFLAGS =
SRC = src/main.c src/output.c src/solver.c
OBJ = $(SRC:.c=.o)
TARGET = bin/cfd_solver.o

all: $(TARGET)

$(TARGET): $(OBJ)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

src/%.o: src/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf src/*.o $(TARGET)