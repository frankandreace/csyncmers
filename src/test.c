#include "csyncmers.h"
#include "benchmarking.h"
#include <time.h>

#define MINIMIZER_VECTOR 50000

int main(int argc, char *argv[]) {

    if(argc <4) {
        fprintf(stderr, "Usage: %s SEQUENCE K S\n", argv[0]);
        return 1;
    }

    char *sequence_input = argv[1];
    int K = atoi(argv[2]);
    int S = atoi(argv[3]);

    if(S >= K) {
        fprintf(stderr, "Error: S must be less than K\n");
        return 1;
    }

    // Start timing

    printf("RUNNING DEQUE\n") ;
    clock_t start_time_deque = clock();
    int num_results_deque;
    MinimizerResult results_deque[MINIMIZER_VECTOR];
    char* sequence_input_deque = sequence_input;
    compute_closed_syncmers_deque(sequence_input_deque, strlen(sequence_input_deque), K, S, results_deque, &num_results_deque);
    clock_t end_time_deque = clock();
    printf("Number of closed syncmers is [DEQUE]: %d\n", num_results_deque);
    printf("END DEQUE\n\n") ;

    printf("RUNNING NAIVE\n") ;
    clock_t start_time_naive = clock();
    int num_results_naive;
    MinimizerResult results_naive[MINIMIZER_VECTOR];
    char* sequence_input_naive = sequence_input;
    compute_closed_syncmers_naive(sequence_input_naive, strlen(sequence_input_naive), K, S, results_naive, &num_results_naive);
    clock_t end_time_naive = clock();
    printf("Number of closed syncmers is [NAIVE]: %d\n", num_results_naive);
    printf("END NAIVE\n\n") ;

    printf("RUNNING SYNG\n") ;
    clock_t start_time_syng = clock();
    char* sequence_input_syng = sequence_input;
    int num_results_syng;
    MinimizerResult results_syng[MINIMIZER_VECTOR];
    compute_closed_syncmers_syng(sequence_input_syng, strlen(sequence_input_syng), K, S, results_syng, &num_results_syng);
    clock_t end_time_syng = clock();
    printf("END SYNG\n\n") ;

   
    print_benchmark("DEQUE",start_time_deque, end_time_deque, num_results_deque, NULL) ;

    print_benchmark("NAIVE",start_time_naive, end_time_naive, num_results_naive, NULL ) ;

    print_benchmark("SYNG",start_time_syng, end_time_syng, num_results_syng, NULL) ;

    return 0;
}

