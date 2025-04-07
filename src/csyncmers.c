/*
*  MODIFIED / RESTRUCTURED VERSION OF RICHARD DURBIN'S CODE IN SYNG
*  AND RAYAN CHIKHI'S CODE IN CSYNCMERS
*-------------------------------------------------------------------
*  SYNG: https://github.com/richarddurbin/syng/
*  Author: Richard Durbin (rd109@cam.ac.uk)
*  Copyright (C) Richard Durbin, Cambridge University, 2018
*-------------------------------------------------------------------
*  CSYNCMERS: https://github.com/rchikhi/csyncmers
*  Author: RAYAN CHIKHI 
*-------------------------------------------------------------------
*/

#include "csyncmers.h"

void compute_closed_syncmers_syng(char *sequence_input, int len, int K, int S) {
    //, MinimizerResult *results, int *num_results
  /* Durbin's SYNG implementation of closed syncmers enrumeration from a sequence*/
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }
    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    // *num_results = 0;

    // Durbin's function works with the parameters as follows:
    // its K is normal S
    // its W is the window size, i.e. K-S+1
    size_t window_size = (U64)K - (U64)S + 1;

    Seqhash *sh = seqhashCreate(S, window_size, seed);

    // initializing the syncmer iterator
    SeqhashIterator *si = syncmerIterator(sh, sequence_input, len); 
    size_t count = 0 ;

    if (si->iMin != U64MAX) {
        U64 kmer ;
        size_t k_pos = U64MAX;
        size_t last_k_pos;
        size_t s_pos ;
        bool isF ;

        while (syncmerNext(si, &kmer, &k_pos, &s_pos, &isF))
        {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            last_k_pos = k_pos ;
            count++;
        }
        // printf("I SEE k_pos == %lu, last_k_pos == %lu\n", k_pos, last_k_pos) ;
        if(k_pos != U64MAX && last_k_pos != k_pos) {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            count++;
        }
    }
    printf("COUNT IS %lu\n", count) ;
    // releasing memory
    seqhashIteratorDestroy(si);
    seqhashDestroy(sh);

    #ifdef DEBUG
    for(int i = 0; i < count; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
    }
    #endif
}


void compute_closed_syncmers_syng_original(char *sequence_input, int len, int K, int S) {
    //, MinimizerResult *results, int *num_results
  /* Durbin's SYNG implementation of closed syncmers enrumeration from a sequence*/
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }
    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    // *num_results = 0;

    // Durbin's function works with the parameters as follows:
    // its K is normal S
    // its W is the window size, i.e. K-S+1
    size_t window_size = (U64)K - (U64)S + 1;

    Seqhash *sh = seqhashCreate(S, window_size, seed);

    // initializing the syncmer iterator
    SeqhashIterator *si = syncmerIterator_original(sh, sequence_input, len); 
    size_t count = 0 ;

    if (si->iMin != U64MAX) {
        U64 kmer ;
        size_t s_pos = U64MAX;
        size_t last_k_pos;
        bool isF ;

        while (syncmerNext_original(si, &kmer, &s_pos, &isF))
        {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            last_k_pos = s_pos ;
            count++;
        }
        // printf("I SEE k_pos == %lu, last_k_pos == %lu\n", k_pos, last_k_pos) ;
        if(s_pos != U64MAX && last_k_pos != s_pos) {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            count++;
        }
    }
    printf("COUNT IS %lu\n", count) ;
    // releasing memory
    seqhashIteratorDestroy(si);
    seqhashDestroy(sh);

    #ifdef DEBUG
    for(int i = 0; i < count; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
    }
    #endif
}

// , MinimizerResult *results, int *num_results
void compute_closed_syncmers_deque_iterator(char *sequence_input, int len, int K, int S){
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }
    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    // *num_results = 0;
    size_t window_size = K - S + 1;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    // printf("SYNCMER DEQUE INITIALIZATION\n") ;
    SeqhashIterator *si = syncmerDequeIterator(sh, sequence_input, len);
    // printf("END SYNCMER DEQUE INITIALIZATION\n") ;
    U64 smer ;
    size_t count = 0 ;
    size_t s_pos ;
    bool isF ;
    // printf("CHECKIN IMIN\n") ;
    // printf("IMIN IS %lu", si->iMin);
    if (si->iMin != U64MAX) {
        U64 kmer ;
        size_t k_pos = U64MAX;
        size_t last_k_pos;
        size_t s_pos ;
        bool isF ;
        // printf("STARTING WHILE\n") ;
        while (syncmerDequeNext(si, &kmer, &k_pos, &s_pos, &isF))
        {
            // printf("ADDING\n") ;
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            last_k_pos = k_pos ;
            count++;
        }
        // printf("I SEE k_pos == %lu, last_k_pos == %lu\n", k_pos, last_k_pos) ;
        if(k_pos != U64MAX && last_k_pos != k_pos) {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            count++;
        }
    }
    // printf("COUNT IS %lu\n", count) ;
    // releasing memory
    seqhashIteratorDestroy(si);
    seqhashDestroy(sh);

}

void compute_closed_syncmers_deque(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results) {
  /* Rayan's implementation of the sliding window approach */
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    *num_results = 0;

    size_t num_s_mers = len - S + 1;
    U64 *s_mer_hashes = (U64 *)malloc(num_s_mers * sizeof(U64));
    size_t window_size = K - S + 1;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, len);

    U64 smer;
    size_t s_pos;
    bool isF;
    int ii = 0;

    // PRE-ENTIVELY HASH ALL S-MERS
    while(seqhashNext(si, &smer, &s_pos, &isF)){
        s_mer_hashes[ii++] = smer;
    }

    // INITIALIZE DEQUE
    size_t *deque = (size_t *)malloc(num_s_mers * sizeof(size_t));
    size_t front = 0, back = 0;
    size_t min_pos, kmer_pos;

    // USE DEQUE TO FIND MINIMAL S-MER IN O(N)
    for(size_t i = 0; i < num_s_mers; i++){
        while(back > front && s_mer_hashes[deque[back-1]] > s_mer_hashes[i]){
            back--;
        }
        deque[back++] = i;
        if(i >= window_size && deque[front] <= i - window_size) {
            front++;
        }

        // CHECKING FOR CLOSED SYNCMER CONDITION
        if (i >= window_size - 1) {
            min_pos = deque[front];
            kmer_pos = i - window_size + 1;

            #ifdef DEBUG
            printf("ITERATION %lu: min_pos: %lu; min_val: %llu; start_pos: %lu, start_val: %llu; end_pos: %lu, end_val: %llu\n", kmer_pos, min_pos, s_mer_hashes[min_pos], kmer_pos, s_mer_hashes[kmer_pos], kmer_pos + K - S, s_mer_hashes[kmer_pos + K - S]) ;
            #endif

            if (min_pos == kmer_pos || min_pos == kmer_pos + K - S){
                add_minimizer(results, num_results, s_mer_hashes[min_pos], kmer_pos, min_pos);
            }

        }

    }

    #ifdef DEBUG
    for(int i = 0; i < num_s_mers; i++){
        printf("%llu\n", s_mer_hashes[i]);
    }
    for(int i = 0; i < *num_results; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
    }
    #endif

    //releasing memory
    free(s_mer_hashes);
    free(deque);
}

void compute_closed_syncmer_naive_iterator(char *sequence_input, int len, int K, int S) {

    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    // setting the seed to 7 as in Durbin's
    U64 seed  = 7 ;
    // *num_results = 0 ;
    size_t length = (size_t)len ;

    // size_t num_s_mers = length - (size_t)S + 1;
    // size_t num_k_mers = length - (size_t)K + 1;
    size_t window_size = (size_t)K - (size_t)S + 1 ;
    
    Seqhash *sh = seqhashCreate(S, window_size, seed) ;

    SeqhashIterator *si = syncmerNaiveIterator(sh, sequence_input, len) ;
    
    U64 smer ;
    size_t count = 0 ;
    size_t s_pos ;
    bool isF ;
    // printf("CHECKIN IMIN\n") ;
    // printf("IMIN IS %lu", si->iMin);
    if (si->iMin != U64MAX) {
        U64 kmer ;
        size_t k_pos = U64MAX;
        size_t last_k_pos;
        size_t s_pos ;
        bool isF ;
        // printf("STARTING WHILE\n") ;
        while (syncmerNaiveNext(si, &kmer, &k_pos, &s_pos, &isF))
        {
            // printf("ADDING\n") ;
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            last_k_pos = k_pos ;
            count++;
        }
        // printf("I SEE k_pos == %lu, last_k_pos == %lu\n", k_pos, last_k_pos) ;
        if(k_pos != U64MAX && last_k_pos != k_pos) {
            // add_minimizer(results, num_results, kmer, k_pos, s_pos);
            count++;
        }
    }
    printf("COUNT IS %lu\n", count) ;
    // releasing memory
    seqhashIteratorDestroy(si);
    seqhashDestroy(sh);

}

void compute_closed_syncmers_naive(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results) {
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    *num_results = 0;
    size_t length = (size_t)len;

    size_t num_s_mers = length - (size_t)S + 1;
    size_t num_k_mers = length - (size_t)K + 1;
    U64 *s_mer_hashes = (U64 *)malloc(num_s_mers * sizeof(U64));
    size_t window_size = (size_t)K - (size_t)S + 1;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, len);

    U64 smer;
    size_t s_pos;
    bool isF;

    // HASH ALL S-MERS
    while(seqhashNext(si, &smer, &s_pos, &isF)){
        printf("HASH COMPUTED IS %llu\n", smer) ;
        s_mer_hashes[s_pos] = smer;
    }

    size_t front = 0, back = 0;
    size_t min_pos;
    U64 min_smer;

    // USE ARRAY SCAN TO COMPUTE SYNCMERS IN O(N*K)

    for(size_t i = 0; i < num_k_mers; i++){
        min_smer = U64MAX;
        for (size_t j = i; j < i + window_size; j++){

            if (s_mer_hashes[j] < min_smer){
                min_pos = j;
                min_smer = s_mer_hashes[min_pos];
            }

        }

        #ifdef DEBUG
        printf("ITERATION %lu: min_pos: %lu; min_val: %llu; start_pos: %lu, start_val: %llu; end_pos: %lu, end_val: %llu\n", i, min_pos, min_smer, i, s_mer_hashes[i], i+K, s_mer_hashes[i+window_size]) ;
        #endif
        
        if (min_pos == i || min_pos == i + window_size-1){
            add_minimizer(results,num_results,min_smer,i,min_pos);
        } 
    }

    // DEBUG PRINTS
    #ifdef DEBUG
    for(int i = 0; i < num_s_mers; i++){
        printf("%llu\n", s_mer_hashes[i]);
    }
    for(int i = 0; i < *num_results; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
    }
    #endif

    // releasing memory
    free(s_mer_hashes) ;
}