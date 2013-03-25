NUM?=256
TESTS=$(shell ls *.cpp | sed -e 's/.cpp$$//g';)

all: $(TESTS)

run: all
	@for b in $(TESTS); do ./$$b; done

../driver.o: ../driver.cpp
	make -C ../ driver.o

clean:
	rm -f $(TESTS) *.o *.s

-include ../def.mk
