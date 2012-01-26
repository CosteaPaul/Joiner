all: join

join:
	g++ join.cpp -fopenmp -O3 -o join

clean:
	rm -rf *o join