/***************************************************************************
 *
 * Parallel version of Gaussian elimination
 *
 ***************************************************************************/

#define _XOPEN_SOURCE 600

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>

#define MAX_SIZE 4096

typedef double matrix[MAX_SIZE][MAX_SIZE];

int N;				/* matrix size		*/
int maxnum;			/* max number of element*/
char *Init;			/* matrix init type	*/
int PRINT;			/* print switch		*/
matrix A;			/* matrix A		*/
double b[MAX_SIZE]; /* vector b             */
double y[MAX_SIZE]; /* vector y             */

int T;				/* thread count		*/
pthread_barrier_t barrier;

/* forward declarations */
void *thread_worker(void *);
void work(void);
void Init_Matrix(void);
void Print_Matrix(void);
void Init_Default(void);
int Read_Options(int, char **);

int main(int argc, char **argv)
{
	int i, timestart, timeend, iter;
    long sec_elapsed, msec_elapsed;
	struct timeval start, stop;

	Init_Default();			  /* Init default values	*/
	Read_Options(argc, argv); /* Read arguments	*/
	Init_Matrix();			  /* Init the matrix	*/

    gettimeofday(&start, NULL);
	work();
    gettimeofday(&stop, NULL);

    sec_elapsed = (stop.tv_sec - start.tv_sec);
    msec_elapsed = ((sec_elapsed*1000000) + stop.tv_usec) - (start.tv_usec);

	if (PRINT == 1)
		Print_Matrix();

	printf("Time taken: %ld ms\n", msec_elapsed / 1000);
}

void *thread_worker(void *arg)
{
	int tid = (long long)arg;
	int i, j, k;

	for (k = 0; k < N; k++)
	{ /* Outer loop */
		for (j = k + 1 + tid; j < N; j += T)
			A[k][j] = A[k][j] / A[k][k]; /* Division step */

		if (tid == 0)
			y[k] = b[k] / A[k][k];

		pthread_barrier_wait(&barrier);

		if (tid == 0)
			A[k][k] = 1.0;

		for (i = k + 1 + tid; i < N; i += T)
		{
			for (j = k + 1; j < N; j++)
				A[i][j] = A[i][j] - A[i][k] * A[k][j]; /* Elimination step */

			b[i] = b[i] - A[i][k] * y[k];
			A[i][k] = 0.0;
		}

		pthread_barrier_wait(&barrier);
	}

    return NULL;
}

void work(void)
{
	printf("threads   = %d \n\n", T);

	pthread_t *threads = malloc(sizeof(pthread_t *) * (T - 1));

	// Initialize barrier
	pthread_barrier_init(&barrier, NULL, T);

	// Create threads
	for (int i = 0; i < T - 1; ++i)
	{
		pthread_create(&threads[i], NULL, thread_worker, (void *)(long long)(i + 1));
	}

	thread_worker((void *)0ll);

	// Join threads
	for (int i = 0; i < T - 1; ++i)
	{
		pthread_join(threads[i], NULL);
	}

	free(threads);

	// Destroy barrier
	pthread_barrier_destroy(&barrier);
}

void Init_Matrix()
{
	int i, j;

	printf("\nsize      = %dx%d ", N, N);
	printf("\nmaxnum    = %d \n", maxnum);
	printf("Init	  = %s \n", Init);
	printf("Initializing matrix...");

	if (strcmp(Init, "rand") == 0)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				if (i == j) /* diagonal dominance */
					A[i][j] = (double)(rand() % maxnum) + 5.0;
				else
					A[i][j] = (double)(rand() % maxnum) + 1.0;
			}
		}
	}
	if (strcmp(Init, "fast") == 0)
	{
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				if (i == j) /* diagonal dominance */
					A[i][j] = 5.0;
				else
					A[i][j] = 2.0;
			}
		}
	}

	/* Initialize vectors b and y */
	for (i = 0; i < N; i++)
	{
		b[i] = 2.0;
		y[i] = 1.0;
	}

	printf("done \n\n");
	if (PRINT == 1)
		Print_Matrix();
}

void Print_Matrix()
{
	int i, j;

	printf("Matrix A:\n");
	for (i = 0; i < N; i++)
	{
		printf("[");
		for (j = 0; j < N; j++)
			printf(" %5.2f,", A[i][j]);
		printf("]\n");
	}
	printf("Vector b:\n[");
	for (j = 0; j < N; j++)
		printf(" %5.2f,", b[j]);
	printf("]\n");
	printf("Vector y:\n[");
	for (j = 0; j < N; j++)
		printf(" %5.2f,", y[j]);
	printf("]\n");
	printf("\n\n");
}

void Init_Default()
{
	N = 2048;
	T = 4;
	Init = "rand";
	maxnum = 15.0;
	PRINT = 0;
}

int Read_Options(int argc, char **argv)
{
	char *prog;

	prog = *argv;
	while (++argv, --argc > 0)
	{
		if (**argv == '-')
		{
			switch (*++*argv)
			{
			case 'n':
				--argc;
				N = atoi(*++argv);
				break;
			case 't':
				--argc;
				T = atoi(*++argv);
				break;
			case 'h':
				printf("\nHELP: try sor -u \n\n");
				exit(0);
				break;
			case 'u':
				printf("\nUsage: gaussian [-n problemsize]\n");
				printf("           [-D] show default values \n");
				printf("           [-h] help \n");
				printf("           [-I init_type] fast/rand \n");
				printf("           [-m maxnum] max random no \n");
				printf("           [-P print_switch] 0/1 \n");
				exit(0);
				break;
			case 'D':
				printf("\nDefault:  n         = %d ", N);
				printf("\n          t      	  = %d", T);
				printf("\n          Init      = rand");
				printf("\n          maxnum    = 5 ");
				printf("\n          P         = 0 \n\n");
				exit(0);
				break;
			case 'I':
				--argc;
				Init = *++argv;
				break;
			case 'm':
				--argc;
				maxnum = atoi(*++argv);
				break;
			case 'P':
				--argc;
				PRINT = atoi(*++argv);
				break;
			default:
				printf("%s: ignored option: -%s\n", prog, *argv);
				printf("HELP: try %s -u \n\n", prog);
				break;
			}
		}
	}
}
