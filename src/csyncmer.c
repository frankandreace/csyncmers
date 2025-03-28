// #include "closed_syncmers.h"
// #include "closed_syncmers_naive.h"
#include "closed_syncmers_syng.h"
#include "fasta_reader.h"
#include <time.h>

#define MINIMIZER_VECTOR 50000
// // Get file size in bytes
// off_t get_file_size(const char *filename) {
//     struct stat st;
//     if (stat(filename, &st) == 0) {
//       return st.st_size;
//     }
//     return -1;
//   }

void print_timing_stats(char* name, double start_time, double end_time, int syncmer_count){

  double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  
  // Calculate and print performance metrics
//   double file_size_mb = file_size / (1024.0 * 1024.0);
//   double processing_speed = file_size_mb / elapsed_time;
  
  printf("%s Performance metrics:\n", name);
//   printf("File size: %.2f MB\n", file_size_mb);
  printf("Processing time: %.2f seconds\n", elapsed_time);
//   printf("Processing speed: %.2f MB/sec\n", processing_speed);
  printf("Syncmers found: %d\n", syncmer_count);
  printf("\n") ;
}

int main(int argc, char *argv[]) {

    if(argc <4) {
        fprintf(stderr, "Usage: %s seqfile K S\n", argv[0]);
        return 1;
    }
    // char *sequence_input = argv[1];
    int K = atoi(argv[2]);
    int S = atoi(argv[3]);
    int isFile = atoi(argv[4]) ;
    if(S >= K) {
        fprintf(stderr, "Error: S must be less than K\n");
        return 1;
    }
    char* sequence_input;

    if (isFile){
    FILE *fin = fopen(argv[1], "r");
    if (!fin) {
        fprintf(stderr, "Error: Cannot open file %s\n", argv[1]);
        return 1;
    }
    stream *s = stream_open_fasta(fin);

    sequence_input = read_sequence(s) ;
    }
    else{
      sequence_input = argv[1] ;
    }
    // Start timing
    clock_t start_time_deque = clock();
  
    int num_results_deque;
    MinimizerResult results_deque[MINIMIZER_VECTOR];
    char* sequence_input_deque = sequence_input;
    compute_closed_syncmers_deque(sequence_input_deque, strlen(sequence_input_deque), K, S, results_deque, &num_results_deque);
    clock_t end_time_deque = clock();
    printf("Number of closed syncmers is [DEQUE]: %d\n", num_results_deque);

    clock_t start_time_naive = clock();
    int num_results_naive;
    MinimizerResult results_naive[MINIMIZER_VECTOR];
    char* sequence_input_naive = sequence_input;
    compute_closed_syncmers_naive(sequence_input_naive, strlen(sequence_input_naive), K, S, results_naive, &num_results_naive);
    clock_t end_time_naive = clock();
    printf("Number of closed syncmers is [NAIVE]: %d\n", num_results_naive);


    clock_t start_time_syng = clock();
    char* sequence_input_syng = sequence_input;
    int num_results_syng;
    MinimizerResult results_syng[MINIMIZER_VECTOR];
    compute_closed_syncmers_syng(sequence_input_syng, strlen(sequence_input_syng), K, S, results_syng, &num_results_syng);
    clock_t end_time_syng = clock();

   
    print_timing_stats("DEQUE",start_time_deque, end_time_deque, num_results_deque) ;

    print_timing_stats("NAIVE",start_time_naive, end_time_naive, num_results_naive) ;

    print_timing_stats("SYNG",start_time_syng, end_time_syng, num_results_syng) ;



    // printf("\nClosed Syncmers (Naive):\n");
    // printf("%-20s %-20s\n", "Position", "Minimizer Hash");
    // for (int i = 0; i < syng_results; i++) {
    //     printf("%-20zu %-20llu\n", syng_results[i].kmer_position, (unsigned long long)syng_results[i].minimizer_hash);
    // }

    // Compare the results
    // if (num_results_syng != num_results_naive) {
    //     printf("\nMismatch in number of closed syncmers: %d (syng) vs %d (naive)\n", num_results_syng, num_results_naive);
    // } else {
    //     printf("\nNumber of closed syncmers matches: %d\n", num_results_naive);
    // }
    // Compare each result
    // int mismatch = 0;
    // for (int i = 0; i < num_results; i++) {
    //     if (results[i].kmer_position != naive_results[i].kmer_position || results[i].minimizer_hash != naive_results[i].minimizer_hash) {
    //         printf("Mismatch at index %d:\n", i);
    //         printf("  Original -> Position: %zu, Hash: %llu\n", results[i].kmer_position, (unsigned long long)results[i].minimizer_hash);
    //         printf("  Naive    -> Position: %zu, Hash: %llu\n", naive_results[i].kmer_position, (unsigned long long)naive_results[i].minimizer_hash);
    //         mismatch = 1;
	//     exit(1);
    //     }
    // }
    // if (!mismatch) {
    //     printf("All closed syncmers match between original and naive method.\n");
    // }

    return 0;
}

