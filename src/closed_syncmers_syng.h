/*
*  MODIFIED / RESTRUCTURED VERSION OF RICHARD DURBIN'S CODE IN SYNG
*  SYNG: https://github.com/richarddurbin/syng/
*  Author: Richard Durbin (rd109@cam.ac.uk)
*  Copyright (C) Richard Durbin, Cambridge University, 2018
*-------------------------------------------------------------------
*/

 #include "utils.h"

 typedef struct {
   int seed ;			/* seed */
   int k ;			/* kmer */
   int w ;			/* window */
   U64 mask ;			/* 2*k bits */
   int shift1, shift2 ;
   U64 factor1, factor2 ;
   U64 patternRC[4] ;		/* one per base */
 } Seqhash ;
 
 typedef struct {
   Seqhash *sh ;
   char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
   U64 h, hRC ;			/* current k-mer values */
   U64 *hash ;			/* buffer of length w holding hashes for current window */
   bool *isForward ;		/* buffer of length w holding isForward for current window */
   int base ;			/* start of buf in sequence */
   size_t iStart, iMin ;		/* position in buf of start of current window, next min */
   U64 min ;                     /* needed for syncmers */
   bool isDone ;
 } SeqhashIterator ;
 
 Seqhash *seqhashCreate (int k, int w, int seed) ;
 static void seqhashDestroy (Seqhash *sh) { free (sh) ; }
 
// simple iterator to return all kmers
SeqhashIterator *seqhashIterator (Seqhash *sh, char *s, int len) ;
bool seqhashNext (SeqhashIterator *si, U64 *kmer, size_t *pos, bool *isF) ;

static void seqhashIteratorDestroy (SeqhashIterator *si)
{ free (si->hash) ; free (si->isForward) ; free (si) ; }

static inline U64 hashRC (SeqhashIterator *si, bool *isForward);
static inline U64 advanceHashRC (SeqhashIterator *si, bool *isForward);
// (closed) syncmer extracts w-mers that end with a minimal kmer
// these provide a cover, and have good distribution properties
SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len) ;
bool syncmerNext (SeqhashIterator *si, U64 *kmer, size_t *k_pos, bool *isF) ;

// utilities
static inline U64 kHash (Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1) ; }

// wrap function
void compute_closed_syncmers_syng(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

void compute_closed_syncmers_deque(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

void compute_closed_syncmers_naive(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results);

/******* end of file ********/