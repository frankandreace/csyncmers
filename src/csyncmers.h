#include "hashing.h"
#include "minimizer_storage.h"


// wrap functions to compute closed syncmers
void compute_closed_syncmers_syng(char *sequence_input, int len, int K, int S);

void compute_closed_syncmers_deque(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

void compute_closed_syncmers_deque_iterator(char *sequence_input, int len, int K, int S);
// , MinimizerResult *results, int *num_results
void compute_closed_syncmers_naive(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

void compute_closed_syncmers_syng_original(char *sequence_input, int len, int K, int S);

void compute_closed_syncmer_naive_iterator(char *sequence_input, int len, int K, int S);
/******* end of file ********/