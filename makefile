#replace the .dylibs with the files you want to make, and replace /usr/local/bin/gcc-5 with your version of gcc. 


all:
	$(MAKE) -C usefulfunctions
	$(MAKE) -C rk4
	$(MAKE) -C sparse_matrix_table
	$(MAKE) -C invariantdynamics 
	$(MAKE) -C  figurefactory
clean:
	rm -f */*.so */*.o */*.dylib */*/*.dylib

