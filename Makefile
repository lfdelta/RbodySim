
tests:
	make rbody
	make sim
	make continuous

rbody:
	clang++ testrbody.cpp -Wall -o testrbody

sim:
	clang++ testsim.cpp -Wall -Wno-c++11-extensions -Wno-char-subscripts -o testsim

continuous:
	clang++ testcontinuous.cpp -Wall -Wno-c++11-extensions -Wno-char-subscripts -o testcont

clean:
	rm ./testsim ; rm ./testrbody ; rm ./testcont ; rm *~ *.dyn *.sys *.sta *.sim