all: rk4.dylib
exe: rk4

%.dylib : %.c
	/usr/local/bin/gcc-5 -O3 -std=c99 -lm -lquadmath -dynamiclib -o $@ $^


#% : %.c 
#	/usr/local/bin/gcc-5 -O3 -std=c99 -o $@ $^

#%.o : %.c
#	/usr/local/bin/gcc-5 -O3 -std=c99 -fPIC -c $^

#%.dylib : %.o
#	/usr/local/bin/gcc-5 -dynamiclib -o $@ $^

#%.so: %.o
#	/usr/local/bin/gcc-5 -shared -o $@ $^

clean:
	rm -f *.so *.o *.dylib

