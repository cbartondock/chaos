all: ergodic.dylib  invariant4.dylib quasiperiodicity.dylib 
exe: ergodic  invariant4 quasiperiodicity

% : %.c 
	/usr/local/bin/gcc-5 -O3 -std=c99 -fopenmp -o $@ $^

%.o : %.c
	/usr/local/bin/gcc-5 -O3 -std=c99 -fopenmp -fPIC -c $^

%.dylib : %.o
	/usr/local/bin/gcc-5 -dynamiclib -fopenmp -o $@ $^

%.so: %.o
	/usr/local/bin/gcc-5 -shared -o $@ $^

clean:
	rm -f *.so *.o *.dylib

