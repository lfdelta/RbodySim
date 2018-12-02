
tests:
	make rbody
	make euler
	make continuous

rbody:
	clang++ testrbody.cpp -o testrbody

euler:
	clang++ testeuler.cpp -o testeuler

continuous:
	clang++ testcontinuous.cpp -o testcont

clean:
	rm ./testeuler ; rm ./testrbody ; rm *~ *.dyn *.sys *.sta *.sim