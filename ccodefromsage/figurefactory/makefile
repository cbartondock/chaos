all: test.dylib
exe: test


% : %.c 
	/opt/local/bin/gcc-mp-4.9 -O3 -std=c99 -fopenmp -o $@ $^

%.o : %.c
	/opt/local/bin/gcc-mp-4.9 -O3 -std=c99 -fopenmp -fPIC -c $^

%.dylib : %.o
	/opt/local/bin/gcc-mp-4.9 -dynamiclib -fopenmp -o $@ $^

%.so: %.o
	/opt/local/bin/gcc-mp-4.9 -shared -o $@ $^

clean:
	rm -f *.so *.o *.dylib

