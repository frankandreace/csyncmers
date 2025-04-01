#include "hashing.h"
#include "minimizer_storage.h"


// wrap functions to compute closed syncmers
void compute_closed_syncmers_syng(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

void compute_closed_syncmers_deque(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

void compute_closed_syncmers_naive(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

/******* end of file ********/