
gauss: gauss_sequential gauss_parallel

qsort: qsort_sequential qsort_parallel

all: gauss_sequential gauss_parallel qsort_sequential qsort_parallel

gauss_sequential: gaussianseq.c
	gcc -O2 -o gauss_sequential gaussianseq.c -lpthread

gauss_parallel: gaussianpar.c
	gcc -O2 -o gauss_parallel gaussianpar.c -lpthread

qsort_sequential: qsortseq.c
	gcc -O2 -o qsort_sequential qsortseq.c -lpthread

qsort_parallel: qsortpar.c
	gcc -O2 -o qsort_parallel qsortpar.c -lpthread

clean:
	rm *.o *.a gauss_sequential gauss_parallel qsort_sequential qsort_parallel
