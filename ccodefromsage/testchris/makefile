all: test.dylib
exe: test


% : %.c 
	gcc -O3 -std=c99 -o $@ $^

%.o : %.c
	gcc -O3 -std=c99 -fPIC -c $^

%.dylib : %.o
	gcc -dynamiclib -o $@ $^

%.so: %.o
	gcc -shared -o $@ $^

clean:
	rm -f *.so *.o *.dylib

