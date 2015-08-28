
FLAGS=-O2 -shared -fPIC

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
	gcc -L./ -Wl,-rpath=./ ChiSqTester.c -O2 -o tester -lChiSq
	@echo ""

libChiSq.so: ChiSquare.c ChiSquare.h
	@echo ""
	@echo "Building shared library"
	@echo ""
	gcc $(FLAGS) ChiSquare.c -o libChiSq.so
	@echo ""

.PHONY:clean
clean:
	-rm libChiSq.so
	-rm tester
