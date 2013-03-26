FOLDERS=$(shell ls -d *-kernels;)
TARGETS=$(shell ls -d *-kernels |  sed -e 's/-kernels$$//g';)

all: $(TARGETS)

run: all
	@for i in $(TARGETS); do echo "STREAM"; make -s -C $$i-kernels run; done

driver.o: driver.cpp
	g++ -DHRT_ARCH=2 -I../../copy_benchmark/hrtimer $(CPPFLAGS) -c -o $@ $< 

%: %-kernels %-kernels/*.cpp driver.o
	@make -C $<

clean:
	rm -f *.o
	@for i in $(FOLDERS); do make -C $$i clean; done

-include def.mk
