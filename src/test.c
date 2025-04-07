#include "csyncmers.h"
#include "benchmarking.h"
#include <time.h>

#define MINIMIZER_VECTOR 50000


int compare_results(MinimizerResult *result_deque, MinimizerResult *result_naive, MinimizerResult *result_syng, int num_deque, int num_naive, int num_syng){
  if (num_deque != num_naive || num_syng != num_deque ){
    printf("ERROR! NUMBER OF RESULT IS NOT CONCORDANT.\n") ;
    printf("NUM DEQUE: %d ; NUM NAIVE: %d ; NUM SYNG : %d .\n", num_deque, num_naive, num_syng) ;
    dump_to_disk( result_deque, num_deque, "./result_deque.txt") ;
    dump_to_disk( result_naive, num_naive, "./result_naive.txt") ;
    dump_to_disk( result_syng, num_syng, "./result_syng.txt") ;
    return 1 ;
  }

  int mismatch = 0 ;

  for( int i = 0 ; i < num_deque; i++){
    if (result_deque[i].minimizer_hash != result_naive[i].minimizer_hash || result_deque[i].minimizer_hash != result_syng[i].minimizer_hash){
      printf("ERROR ON %d-th syncmer minimizer.\n", i) ;
      printf("DEQUE MINIMIZER: %llu; NAIVE MINIMIZER: %llu; SYNG MINIMIZER: %llu\n", result_deque[i].minimizer_hash, result_naive[i].minimizer_hash, result_syng[i].minimizer_hash) ;
      dump_to_disk( result_deque, num_deque, "./result_deque.txt") ;
      dump_to_disk( result_naive, num_naive, "./result_naive.txt") ;
      dump_to_disk( result_syng, num_syng, "./result_syng.txt") ;
      return 1 ;
    }
    if (result_deque[i].kmer_position != result_naive[i].kmer_position || result_deque[i].kmer_position != result_syng[i].kmer_position){
      printf("ERROR ON %d-th syncmer k-mer position.\n", i) ;
      printf("DEQUE K-MER POSITION: %lu; NAIVE K-MER POSITION: %lu; SYNG K-MER POSITION: %lu\n", result_deque[i].kmer_position, result_naive[i].kmer_position, result_syng[i].kmer_position) ;
      dump_to_disk( result_deque, num_deque, "./result_deque.txt") ;
      dump_to_disk( result_naive, num_naive, "./result_naive.txt") ;
      dump_to_disk( result_syng, num_syng, "./result_syng.txt") ;
      return 1 ;
    }
    if (result_deque[i].smer_position != result_naive[i].smer_position || result_deque[i].smer_position != result_syng[i].smer_position){
      printf("ERROR ON %d-th syncmer s-mer position.\n", i) ;
      printf("DEQUE S-MER POSITION: %lu; NAIVE S-MER POSITION: %lu; SYNG S-MER POSITION: %lu\n", result_deque[i].smer_position, result_naive[i].smer_position, result_syng[i].smer_position) ;
      dump_to_disk( result_deque, num_deque, "./result_deque.txt") ;
      dump_to_disk( result_naive, num_naive, "./result_naive.txt") ;
      dump_to_disk( result_syng, num_syng, "./result_syng.txt") ;
      return 1 ;
    }
  }
  return 0 ;
}

void dump_to_disk(const MinimizerResult *result, const int num_elements, const char* filename){
  FILE *fp = fopen(filename, "w");
  if (fp == NULL) {
    perror("Error opening file");
    exit(EXIT_FAILURE);
  }

  for(int i = 0 ; i < num_elements ; i++) {
    fprintf(fp, "%llu\t%lu\t%lu\n", result[i].minimizer_hash, result[i].kmer_position, result[i].smer_position) ;
  }

  fclose(fp) ;
  printf("DUMPED to %s\n", filename);
}

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
    int sequence_input_length = strlen(sequence_input);

    //allocating memory
    printf("ALLOCATING MEMORY\n") ;
    char *sequence_input_deque = malloc(sequence_input_length + 1) ;
    char *sequence_input_naive = malloc(sequence_input_length + 1) ; 
    char *sequence_input_syng = malloc(sequence_input_length + 1) ;

    if (sequence_input_deque == NULL) {
      perror("Failed to allocate memory");
      return 1;
    }
    if (sequence_input_naive == NULL) {
      perror("Failed to allocate memory");
      return 1;
    }
    if (sequence_input_syng == NULL) {
      perror("Failed to allocate memory");
    return 1;
    }

    // Copying the string
    printf("COPYING THE STRING\n") ;
    strcpy(sequence_input_deque, sequence_input) ;
    strcpy(sequence_input_naive, sequence_input) ;
    strcpy(sequence_input_syng, sequence_input) ;
    // Start timing

    printf("RUNNING NAIVE\n") ;
    clock_t start_time_naive = clock();
    int num_results_naive;
    MinimizerResult results_naive[MINIMIZER_VECTOR];
    compute_closed_syncmers_naive(sequence_input_naive, sequence_input_length, K, S, results_naive, &num_results_naive);
    clock_t end_time_naive = clock();
    printf("Number of closed syncmers is [NAIVE]: %d\n", num_results_naive);
    printf("END NAIVE\n\n") ;

    printf("RUNNING DEQUE\n") ;
    clock_t start_time_deque = clock();
    int num_results_deque;
    MinimizerResult results_deque[MINIMIZER_VECTOR];

    compute_closed_syncmers_deque_iterator(sequence_input_deque, sequence_input_length, K, S); //, results_deque, &num_results_deque
    clock_t end_time_deque = clock();
    printf("Number of closed syncmers is [DEQUE]: %d\n", num_results_deque);
    printf("END DEQUE\n\n") ;
    // printf("RUNNING SYNG\n") ;
    // clock_t start_time_syng = clock();
    // int num_results_syng;
    // MinimizerResult results_syng[MINIMIZER_VECTOR];
    // compute_closed_syncmers_syng(sequence_input_syng, sequence_input_length, K, S, results_syng, &num_results_syng);
    // clock_t end_time_syng = clock();
    // printf("Number of closed syncmers is [SYNG]: %d\n", num_results_syng);
    // printf("END SYNG\n\n") ;
    // print_benchmark("DEQUE",start_time_deque, end_time_deque, num_results_deque) ;

    // print_benchmark("NAIVE",start_time_naive, end_time_naive, num_results_naive) ;

    // print_benchmark("SYNG",start_time_syng, end_time_syng, num_results_syng, NULL) ;


    // Compare each result
    // compare_results(results_deque, results_naive, results_syng, num_results_deque, num_results_naive, num_results_syng) ;

    // free sequences

    // free(sequence_input);
    // free(sequence_input_deque);
    // free(sequence_input_naive);
    // free(sequence_input_syng);

    // free minimizers
    // free(results_deque);
    // free(results_naive);
    // free(results_syng);

    return compare_results(results_deque, results_naive, results_naive, num_results_deque, num_results_naive, num_results_naive) ;
}

