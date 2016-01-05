all: smtable.dylib  
exe: smtable smtests



%.dylib : %.c
	/usr/local/bin/gcc-5 -O3 -std=c99 -Wall -lm -lquadmath -dynamiclib -o $@ $^

#% : %.c 
#	/usr/local/bin/gcc-5 -O3 -std=c99 -Wall -o -lm -lquadmath $@ $^

#%.o : %.c
#	/usr/local/bin/gcc-5 -O3 -std=c99 -Wall -fPIC -c -lm -lquadmath $^

#%.dylib : %.o
#	/usr/local/bin/gcc-5 -dynamiclib -Wall -o -lm -lquadmath $@ $^
#
#%.so: %.o
#	/usr/local/bin/gcc-5 -shared -o -Wall -lm -lquadmath $@ $^

clean:
	rm -f *.so *.o *.dylib

