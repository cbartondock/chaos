all: dylibs/birkhoff_partition.dylib dylibs/invariant4.dylib dylibs/quasiperiodicity.dylib dylibs/switch.dylib dylibs/sticky.dylib dylibs/quasi_curve.dylib dylibs/extended_smap.dylib dylibs/split_quasi.dylib dylibs/lyapunov.dylib dylibs/lyapunovconv.dylib dylibs/wconv.dylib dylibs/birk_curve.dylib dylibs/optquasiperiodicity.dylib
exe: birkhoff_partition invariant invariant4 quasiperiodicity switch sticky extended_smap weightsum drift_detect rotation_num




dylibs/%.dylib : %.c
	/usr/local/bin/gcc-5 -O3 -std=c99 -Wall -fopenmp -lm -lquadmath -dynamiclib -o $@ $^

clean:
	rm -f *.so *.o dylibs/*.dylib

