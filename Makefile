
tests:
	make rbody
	make sim
	make continuous

rbody:
	clang++ testrbody.cpp -o testrbody

sim:
	clang++ testsim.cpp -o testsim

continuous:
	clang++ testcontinuous.cpp -o testcont

clean:
	rm ./testsim ; rm ./testrbody ; rm ./testcont ; rm *~ *.dyn *.sys *.sta *.sim