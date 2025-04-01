#include "utils.h"

// Dynamically resize array for minimizers
// Struct to hold minimizer and its position
typedef struct {
    U64 minimizer_hash;
    size_t kmer_position;
    size_t smer_position;
    // bool isForward;
} MinimizerResult;

static void add_minimizer(MinimizerResult *results, int *size, U64 minimizer_hash, size_t kmer_position, size_t smer_position) {
    
    #ifdef DEBUG
    printf("ADDING: %llu from kmer position %lu and smer position %lu.\n", minimizer_hash, kmer_position, smer_position) ;
    #endif

    results[*size].minimizer_hash = minimizer_hash;
    results[*size].kmer_position = kmer_position;
    results[*size].smer_position = smer_position;
    (*size)++;
    
    //add resize case
    
}
