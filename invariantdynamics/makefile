all: dylibs/birkhoff_partition.dylib  dylibs/invariant.dylib dylibs/invariant4.dylib dylibs/quasiperiodicity.dylib dylibs/switch.dylib dylibs/sticky.dylib dylibs/quasi_curve.dylib dylibs/extended_smap.dylib dylibs/split_quasi.dylib dylibs/lyapunov.dylib dylibs/lyapunovconv.dylib dylibs/wconv.dylib dylibs/birk_curve.dylib
exe: birkhoff_partition invariant invariant4 quasiperiodicity switch sticky extended_smap weightsum drift_detect

% : %.c 
	/usr/local/bin/gcc-5 -O3 -Wall -fopenmp -std=c99 -o $@ $^

%.o : %.c
	/usr/local/bin/gcc-5 -O3 -Wall -fopenmp -std=c99 -fPIC -c $^

dylibs/%.dylib : %.o
	/usr/local/bin/gcc-5 -Wall -fopenmp -dynamiclib -o $@ $^

%.so: %.o
	/usr/local/bin/gcc-5 -Wall -fopenmp -shared -o $@ $^

clean:
	rm -f *.so *.o dylibs/*.dylib

