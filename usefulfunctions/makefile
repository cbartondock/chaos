all: functions.dylib 
exe: functions

%.dylib : %.c
	/usr/local/bin/gcc-5 -O3 -std=c99 -lm -lquadmath -dynamiclib -o $@ $^

#% : %.c 
#	/usr/local/bin/gcc-5 -O3 -std=c99 -o -lm -lquadmath $@ $^

#%.o : %.c
#	/usr/local/bin/gcc-5 -O3 -std=c99 -fPIC -c -lm -lquadmath $^

#%.dylib : %.o
#	/usr/local/bin/gcc-5 -dynamiclib -o -lm -lquadmath $@ $^

#%.so: %.o
#	/usr/local/bin/gcc-5 -shared -o -lm -lquadmath $@ $^

clean:
	rm -f *.so *.o *.dylib

