all:
	make -C tests
	make -C pmpi-tests

clean:
	make -C tests clean
	make -C pmpi-tests clean
	svn status
