# starfix: makefile containing code formatting, compilation, and testing targets

CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfloat-equal -Werror -O2
DEPS_CFLAGS = -Wall -O2
LDFLAGS = -lm

.PHONY: all indent check clean demo lint db

all: t/test_ekf t/test_solver t/test_centroid t/test_attitude t/test_graph t/test_identify starfix_cli

nanoqsp.o: deps/nanoqsp/nanoqsp.c
	$(CC) $(DEPS_CFLAGS) -c $< -o $@

nanoqsp_blas.o: deps/nanoqsp/nanoqsp_blas.c
	$(CC) $(DEPS_CFLAGS) -c $< -o $@

t/test_ekf: src/starfix_ekf.c src/test_ekf.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_solver: src/starfix_solver.c src/test_solver.c nanoqsp.o nanoqsp_blas.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_centroid: src/starfix_centroid.c src/test_centroid.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_attitude: src/starfix_attitude.c src/test_attitude.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_graph: src/starfix_graph.c src/test_graph.c nanoqsp.o nanoqsp_blas.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

t/test_identify: src/starfix_identify.c src/starfix_attitude.c src/test_identify.c
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

starfix_cli: src/starfix_cli.c src/starfix_centroid.c src/starfix_identify.c src/starfix_attitude.c src/starfix_solver.c src/starfix_ekf.c src/starfix_graph.c nanoqsp.o nanoqsp_blas.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

indent:
	clang-format -i src/*.c src/*.h

lint:
	cppcheck --enable=all --suppress=missingIncludeSystem src/*.c src/*.h
	clang-tidy src/*.c src/*.h -- -Isrc

check: all
	prove -v t/*.t

demo: starfix_cli
	perl scripts/simulate_camera.pl
	./starfix_cli --pipeline
	perl scripts/visualize.pl

db:
	$(MAKE) -C data hyg
	perl scripts/generate_database.pl

clean:
	rm -f t/test_ekf t/test_solver t/test_centroid t/test_attitude t/test_graph t/test_identify starfix_cli *.o
	rm -f data/step1_centroids.png data/step2_identified.png data/step2b_zenith_camera.png data/step2c_zenith_celestial.png data/step2d_earth_coordinates.png data/step3_position.png data/solved_position.json
