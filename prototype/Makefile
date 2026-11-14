CC = gcc
CFLAGS = -Wall -Wextra -O2 -I../src
LDFLAGS = -lm

HW_SRCS = main.c camera.c bno055.c
HW_OBJS = $(HW_SRCS:.c=.o)
TARGET = starfix_prototype

SIM_SRCS = sim.c camera_mock.c bno055_mock.c
SIM_OBJS = $(SIM_SRCS:.c=.o)
SIM_TARGET = starfix_prototype_sim

LIBS = ../libstarfix.a ../libstarfix_navigation.a

.PHONY: all sim check clean libraries

all: $(TARGET) $(SIM_TARGET)

libraries:
	$(MAKE) -C .. libstarfix.a libstarfix_navigation.a

# Check the root dependency graph once, then refresh archive timestamps without
# touching them so unchanged libraries do not force either executable to relink.
$(LIBS): libraries
	@:

$(TARGET): $(HW_OBJS) $(LIBS)
	$(CC) $(HW_OBJS) $(LIBS) -o $(TARGET) $(LDFLAGS)

$(SIM_TARGET): $(SIM_OBJS) $(LIBS)
	$(CC) $(SIM_OBJS) $(LIBS) -o $(SIM_TARGET) $(LDFLAGS)

sim: $(SIM_TARGET)
	./$(SIM_TARGET)

test_camera: test_camera.o
	$(CC) test_camera.o -o $@ $(LDFLAGS)

check: test_camera
	./test_camera
	prove -v ../t/04-prototype-build.t

%.o: %.c
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

-include *.d

clean:
	rm -f $(HW_OBJS) $(SIM_OBJS) test_camera.o *.d $(TARGET) $(SIM_TARGET) test_camera bodega_to_marquesas.tsv bodega_to_marquesas.nmea
