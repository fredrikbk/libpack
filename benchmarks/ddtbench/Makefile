include Makefile.inc
#some parameter for the test environment, in a different file, so that the shells scripts can red them properly
include test.inc
all: ddtbench_c ddtbench_f90

clean:
	$(MAKE) -C src_f90 clean
	$(MAKE) -C src_c clean
	rm -f ddtbench_c ddtbench_f90

distclean: clean
	$(MAKE) -C src_f90 distclean
	$(MAKE) -C src_c distclean
	rm -rf src_mpi
	rm -rf local

ddtbench_c:
	$(MAKE) -C src_c ddtbench
	mv src_c/ddtbench ddtbench_c

# the fortran version depends on the c version
# because all timing stuff is done in c

ddtbench_f90: ddtbench_c
	$(MAKE) -C src_f90 ddtbench
	mv src_f90/ddtbench ddtbench_f90

######################################################################################################
# part of the makefile for the fullrun test

test: helper_scripts/fullrun.sh build_mpis
	helper_scripts/fullrun.sh $(NPROCS) $(OUTER) $(INNER)

# get major version (because the download url requires it)
OMPI_STABLE_MAJOR=$(shell echo $(OMPI_STABLE) | sed -e 's/\([0-9].[0-9]\).*/\1/';)
OMPI_UNSTABLE_MAJOR=$(shell echo $(OMPI_UNSTABLE) | sed -e 's/\([0-9].[0-9]\).*/\1/';)

build_mpis: openmpi_stable openmpi_unstable mpich2_stable mpich2_unstable

mpiclean:
	rm -rf local/*
	$(MAKE) -C src_mpi/openmpi-$(OMPI_STABLE)/ clean
	$(MAKE) -C src_mpi/openmpi-$(OMPI_UNSTABLE)/ clean
	$(MAKE) -C src_mpi/mpich2-$(MPICH2_STABLE)/ clean
	$(MAKE) -C src_mpi/mpich2-$(MPICH2_UNSTABLE)/ clean

openmpi_unstable: local/openmpi-$(OMPI_UNSTABLE)/bin/mpif90 local/openmpi-$(OMPI_UNSTABLE)/bin/mpiexec

local/openmpi-$(OMPI_UNSTABLE)/bin/mpif90 local/openmpi-$(OMPI_UNSTABLE)/bin/mpiexec: src_mpi/openmpi-$(OMPI_UNSTABLE).tar.bz2
	helper_scripts/install_ompi.sh openmpi-$(OMPI_UNSTABLE)

src_mpi/openmpi-$(OMPI_UNSTABLE).tar.bz2:
	@echo "Trying to download Open MPI version $(OMPI_UNSTABLE) from http://www.open-mpi.org/software/ompi/v$(OMPI_UNSTABLE_MAJOR)/downloads/openmpi-$(OMPI_UNSTABLE).tar.bz2"
	@mkdir -p src_mpi
	@cd src_mpi ; \
	wget -q http://www.open-mpi.org/software/ompi/v$(OMPI_UNSTABLE_MAJOR)/downloads/openmpi-$(OMPI_UNSTABLE).tar.bz2 ; \
	cd ..

openmpi_stable: local/openmpi-$(OMPI_STABLE)/bin/mpif90 local/openmpi-$(OMPI_STABLE)/bin/mpiexec

local/openmpi-$(OMPI_STABLE)/bin/mpif90 local/openmpi-$(OMPI_STABLE)/bin/mpiexec: src_mpi/openmpi-$(OMPI_STABLE).tar.bz2
	helper_scripts/install_ompi.sh openmpi-$(OMPI_STABLE)

src_mpi/openmpi-$(OMPI_STABLE).tar.bz2:
	@echo "Trying to download Open MPI version $(OMPI_STABLE) from http://www.open-mpi.org/software/ompi/v$(OMPI_STABLE_MAJOR)/downloads/openmpi-$(OMPI_STABLE).tar.bz2"
	@mkdir -p src_mpi
	@cd src_mpi ; \
	wget -q http://www.open-mpi.org/software/ompi/v$(OMPI_STABLE_MAJOR)/downloads/openmpi-$(OMPI_STABLE).tar.bz2 ; \
	cd ..

mpich2_stable: local/mpich2-$(MPICH2_STABLE)/bin/mpif90 local/mpich2-$(MPICH2_STABLE)/bin/mpiexec

local/mpich2-$(MPICH2_STABLE)/bin/mpif90 local/mpich2-$(MPICH2_STABLE)/bin/mpiexec: src_mpi/mpich2-$(MPICH2_STABLE).tar.gz
	helper_scripts/install_mpich2.sh mpich2-$(MPICH2_STABLE)

src_mpi/mpich2-$(MPICH2_STABLE).tar.gz:
	@echo "Trying to download MPICH2 version $(MPICH2_STABLE) from http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/$(MPICH2_STABLE)/mpich2-$(MPICH2_STABLE).tar.gz"
	@mkdir -p src_mpi
	@cd src_mpi ; \
	wget -q http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/$(MPICH2_STABLE)/mpich2-$(MPICH2_STABLE).tar.gz ; \
	cd ..

mpich2_unstable: local/mpich2-$(MPICH2_UNSTABLE)/bin/mpif90 local/mpich2-$(MPICH2_UNSTABLE)/bin/mpiexec

local/mpich2-$(MPICH2_UNSTABLE)/bin/mpif90 local/mpich2-$(MPICH2_UNSTABLE)/bin/mpiexec: src_mpi/mpich2-$(MPICH2_UNSTABLE).tar.gz
	helper_scripts/install_mpich2.sh mpich2-$(MPICH2_UNSTABLE)

src_mpi/mpich2-$(MPICH2_UNSTABLE).tar.gz:
	@echo "Trying to download MPICH2 version $(MPICH2_UNSTABLE) from http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/$(MPICH2_UNSTABLE)/mpich2-$(MPICH2_UNSTABLE).tar.gz"
	@mkdir -p src_mpi
	@cd src_mpi ; \
	wget http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/$(MPICH2_UNSTABLE)/mpich2-$(MPICH2_UNSTABLE).tar.gz ; \
	cd ..
