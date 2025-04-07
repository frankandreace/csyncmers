#include "csyncmers.h"
#include "benchmarking.h"
#include "fasta_reader.h"
#include <time.h>


int main(int argc, char *argv[]) {

    if(argc <4) {
        fprintf(stderr, "Usage: %s SEQUENCE K S\n", argv[0]);
        return 1;
    }

    char *fasta_filename = argv[1];
    int K = atoi(argv[2]);
    int S = atoi(argv[3]);

    if(S >= K) {
        fprintf(stderr, "Error: S must be less than K\n");
        return 1;
    }

    FILE *seqFile;
    seqFile = fopen(fasta_filename,"r");
    stream *seqStream = stream_open_fasta(seqFile) ;
    char *seq = read_sequence(seqStream) ;
    int sequence_input_length = strlen(seq) ;
    printf("SEQUENCE LENGTH IS %d\n", sequence_input_length) ;

    clock_t start_time = clock();
    compute_closed_syncmer_naive_iterator(seq, sequence_input_length, K, S);
    clock_t end_time = clock();
    char * method_name = "naive" ;
    print_benchmark(method_name, start_time, end_time, fasta_filename) ;
    // double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // printf("ELAPSED TIME: %f\n", elapsed_time) ;
    fclose(seqFile) ;

    seqFile = fopen(fasta_filename,"r");
    seqStream = stream_open_fasta(seqFile) ;
    seq = read_sequence(seqStream) ;
    sequence_input_length = strlen(seq) ;
    printf("SEQUENCE LENGTH IS %d\n", sequence_input_length) ;

    start_time = clock();
    compute_closed_syncmers_syng(seq, sequence_input_length, K, S);
    end_time = clock();
    method_name = "syng - rescan" ;
    print_benchmark(method_name, start_time, end_time, fasta_filename) ;

    fclose(seqFile) ;
    seqFile = fopen(fasta_filename,"r");
    seqStream = stream_open_fasta(seqFile) ;
    seq = read_sequence(seqStream) ;
    sequence_input_length = strlen(seq) ;
    printf("SEQUENCE LENGTH IS %d\n", sequence_input_length) ;

    start_time = clock();
    compute_closed_syncmers_syng_original(seq, sequence_input_length, K, S);
    end_time = clock();
    method_name = "syng - ORIGINAL" ;
    print_benchmark(method_name, start_time, end_time, fasta_filename) ;

    start_time = clock();
    compute_closed_syncmers_deque_iterator(seq, sequence_input_length, K, S);
    end_time = clock();
    method_name = "DEQUE - ITERATOR" ;
    print_benchmark(method_name, start_time, end_time, fasta_filename) ;
    // elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // printf("ELAPSED TIME: %f\n", elapsed_time) ;
    fclose(seqFile) ;

    return 0;
}