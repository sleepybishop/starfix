CC=gcc
CFLAGS=-Wall -Wextra -O2 -I../src -I../deps/nanoqsp
LDFLAGS=-lm

SRCS=main.c camera.c bno055.c
OBJS=$(SRCS:.c=.o)
TARGET=starfix_prototype

# We link against the core StarFix library files in ../src
CORE_SRCS=../src/starfix_centroid.c \
          ../src/starfix_identify.c \
          ../src/starfix_attitude.c \
          ../src/starfix_solver.c \
          ../src/starfix_ekf.c \
          ../src/starfix_graph.c
CORE_OBJS=$(CORE_SRCS:.c=.o)

# NanoQSP library
NANO_SRCS=../deps/nanoqsp/nanoqsp.c ../deps/nanoqsp/nanoqsp_blas.c
NANO_OBJS=$(NANO_SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS) $(CORE_OBJS) $(NANO_OBJS)
	$(CC) $(OBJS) $(CORE_OBJS) $(NANO_OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(CORE_OBJS) $(NANO_OBJS) $(TARGET)

preprocessor:
	gcc -Wall -O2 trail_preprocessor.c -o trail_preprocessor -lm
