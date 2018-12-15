CC = mpicc
CFLAGS = -O2 -std=c99
LIBS = -lm -Wall -Wextra
OBJDIR = ./bin
HDRDIR = ./headers
SRCDIR = ./src

_OBJ =  median-funcs.o mpi-helper.o main.o vpKnn.o
OBJ = $(patsubst %, $(OBJDIR)/%, $(_OBJ))

_DEPS = median-funcs.h mpi-helper.h vpKnn.h
DEPS = $(patsubst %, $(HDRDIR)/%, $(_DEPS))


mainProgram: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $<  $(CFLAGS)

clean:
	rm -rf ./*.csv