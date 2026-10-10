# starfix: makefile containing code formatting, compilation, and testing targets

CC = gcc
AR = ar
CPPFLAGS = -Isrc -Itools
CFLAGS = -std=c99 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfloat-equal -Werror -O2
DEPS_CFLAGS = -Wall -O2
LDFLAGS = -lm

.PHONY: all tests simulations navigation indent check clean demo lint db

CORE_SRCS = src/starfix_photo.c src/starfix_centroid.c src/starfix_identify.c src/starfix_attitude.c src/starfix_solver.c src/starfix_astrometry.c
CORE_OBJS = $(CORE_SRCS:.c=.o)
NAV_OBJS = src/starfix_ekf.o src/starfix_graph.o src/starfix_nmea.o nanoqsp.o nanoqsp_blas.o
TESTS = t/test_cli_io t/test_photo t/test_ekf t/test_solver t/test_centroid t/test_attitude t/test_graph t/test_identify

all: starfix_cli

tests: all $(TESTS)

simulations: examples/starfix_sim examples/starfix_navigation

navigation: libstarfix_navigation.a examples/starfix_navigation

libstarfix.a: $(CORE_OBJS)
	$(AR) rcs $@ $^

libstarfix_navigation.a: $(NAV_OBJS)
	$(AR) rcs $@ $^

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

-include src/*.d tools/*.d examples/*.d

nanoqsp.o: deps/nanoqsp/nanoqsp.c
	$(CC) $(DEPS_CFLAGS) -c $< -o $@

nanoqsp_blas.o: deps/nanoqsp/nanoqsp_blas.c
	$(CC) $(DEPS_CFLAGS) -c $< -o $@

t/test_cli_io: src/test_cli_io.o tools/starfix_cli_io.o
	$(CC) $^ -o $@ $(LDFLAGS)

t/test_photo: src/test_photo.o libstarfix.a
	$(CC) $^ -o $@ $(LDFLAGS)

t/test_ekf: src/starfix_ekf.c src/test_ekf.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_solver: src/starfix_solver.c src/test_solver.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_centroid: src/starfix_centroid.c src/test_centroid.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_attitude: src/starfix_attitude.c src/starfix_astrometry.c src/test_attitude.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_graph: src/starfix_graph.c src/test_graph.c nanoqsp.o nanoqsp_blas.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_identify: src/starfix_identify.c src/starfix_attitude.c src/starfix_astrometry.c src/test_identify.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

starfix_cli: src/starfix_cli.o tools/starfix_cli_io.o libstarfix.a
	$(CC) $^ -o $@ $(LDFLAGS)

examples/starfix_sim: examples/starfix_sim.o tools/starfix_cli_io.o libstarfix.a
	$(CC) $^ -o $@ $(LDFLAGS)

examples/starfix_navigation: examples/starfix_navigation.o libstarfix_navigation.a
	$(CC) $^ -o $@ $(LDFLAGS)

indent:
	clang-format -i src/*.c src/*.h tools/*.c tools/*.h examples/*.c

lint:
	cppcheck --enable=all --suppress=missingIncludeSystem src/*.c src/*.h
	clang-tidy src/*.c src/*.h -- -Isrc

check: tests simulations
	prove -v t/*.t
	PYTHONDONTWRITEBYTECODE=1 python3 t/test_photo.py

demo: starfix_cli
	$(MAKE) -C demo_photo

db:
	$(MAKE) -C data hyg
	perl scripts/generate_database.pl

clean:
	rm -f t/test_cli_io t/test_photo t/test_ekf t/test_solver t/test_centroid t/test_attitude t/test_graph t/test_identify starfix_cli examples/starfix_sim examples/starfix_navigation *.o *.a
	rm -f src/*.o src/*.d tools/*.o tools/*.d examples/*.o examples/*.d
	rm -f data/step1_centroids.png data/step2_identified.png data/step2b_zenith_camera.png data/step2c_zenith_celestial.png data/step2d_earth_coordinates.png data/step3_position.png data/solved_position.json
