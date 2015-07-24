all: functions.dylib 
exe: functions

% : %.c 
	/usr/local/bin/gcc-5 -O3 -std=c99 -o $@ $^

%.o : %.c
	/usr/local/bin/gcc-5 -O3 -std=c99 -fPIC -c $^

%.dylib : %.o
	/usr/local/bin/gcc-5 -dynamiclib -o $@ $^

%.so: %.o
	/usr/local/bin/gcc-5 -shared -o $@ $^

clean:
	rm -f *.so *.o *.dylib

