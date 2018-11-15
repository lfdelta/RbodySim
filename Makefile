
tests:
	make rbody
	make euler

rbody:
	clang++ testrbody.cpp -o testrbody

euler:
	clang++ testeuler.cpp -o testeuler

clean:
	rm ./test*
	rm *~