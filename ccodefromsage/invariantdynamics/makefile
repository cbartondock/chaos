all: ergodic.dylib  invariant4.dylib 
exe: ergodic  invariant4

% : %.c 
	/usr/bin/gcc -O3 -std=c99 -o $@ $^

%.o : %.c
	/usr/bin/gcc -O3 -std=c99 -fPIC -c $^

%.dylib : %.o
	/usr/bin/gcc -dynamiclib -o $@ $^

%.so: %.o
	/usr/bin/gcc -shared -o $@ $^

clean:
	rm -f *.so *.o *.dylib

