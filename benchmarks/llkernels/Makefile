CPPFLAGS=-DHRT_ARCH=2 -I../../copy_benchmark/hrtimer
.PRECIOUS: %.s %.o

BINARIES=$(shell ls *.ll | grep -v "old" | sed -e 's/.ll$$//g';)

all: $(BINARIES)

run: all
	for b in $(BINARIES); do ./$$b; done

%: %.o %_old.o pack.o
	clang -o $@ $^ 

%.o: %.cpp
	clang $(CPPFLAGS) -c -o $@ $< 

%.o: %.s
	clang -O0 -c -o $@ $< 

%.s: %.ll
	llc -O0 -o $@ $< 

clean:
	rm -f $(BINARIES) *.o *.s

