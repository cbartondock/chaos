all: dylibs/birkhoff_partition.dylib  dylibs/invariant.dylib dylibs/invariant4.dylib dylibs/quasiperiodicity.dylib dylibs/switch.dylib dylibs/sticky.dylib dylibs/quasi_curve.dylib 
exe: birkhoff_partition invariant invariant4 quasiperiodicity switch sticky

% : %.c 
	/usr/local/bin/gcc-5 -O3 -Wall -std=c99 -o $@ $^

%.o : %.c
	/usr/local/bin/gcc-5 -O3 -Wall -std=c99 -fPIC -c $^

dylibs/%.dylib : %.o
	/usr/local/bin/gcc-5 -Wall -dynamiclib -o $@ $^

%.so: %.o
	/usr/local/bin/gcc-5 -Wall -shared -o $@ $^

clean:
	rm -f *.so *.o dylibs/*.dylib

