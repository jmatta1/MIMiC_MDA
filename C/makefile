
FLAGS=-O3 -fstrict-aliasing -Wall -std=gnu99 -mtune=native -march=native
LIB_FLAGS= -ftree-vectorize -shared -fPIC #-Bsymbolic-functions -Wl,-z,now
EXEC_FLAGS=-L./ -Wl,-rpath=./ #,-z,now


.PHONY: all
all: libChiSq.so tester

.PHONY: check
check: tester
	@echo ""
	@echo "Commencing testing of shared library"
	@echo ""
	@echo "Running C based tester"
	@./tester
	@echo ""
	@echo "Running Python based tester"
	@python so_py_tester.py
	@echo ""

tester: ChiSqTester.c libChiSq.so ChiSquare.c ChiSquare.h
	@echo ""
	@echo "Building tester for shared library"
	@echo ""
	gcc $(EXEC_FLAGS) $(FLAGS) ChiSqTester.c -o tester -lChiSq -lrt
	@echo ""

libChiSq.so: ChiSquare.c ChiSquare.h
	@echo ""
	@echo "Building shared library"
	@echo ""
	gcc $(LIB_FLAGS) $(FLAGS) ChiSquare.c -o libChiSq.so
	@echo ""

.PHONY:clean
clean:
	-rm libChiSq.so
	-rm tester
