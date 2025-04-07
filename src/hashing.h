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
    U64 *hashVector ;
    size_t current_position;
    size_t window_size;
  } CircularArray ;

  typedef struct {
    U64 *hashVector ;
    size_t *idVector ;
    size_t front ;
    size_t back ;
    size_t window_size ;
    size_t used_pos ;
    size_t smer_pos ;
  } Deque ;
  
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
    size_t abs_k_pos;
    size_t counter;
    CircularArray *ca;
    Deque *dq;
  } SeqhashIterator ;

// typedef struct {
//   Seqhash *sh ;
//   char *s, *sEnd ;
//   size_t iStart, iMin ;
//   U64 min ;                     
//   bool isDone ;
//   size_t abs_k_pos;
//   CircularArray *ca;
// } myIterator ;


// circular array to handle inserted values
CircularArray *circularArrayCreate(size_t size) ;
static void circularArrayDestroy (CircularArray *ca) { free (ca) ; }
void circularInsert(CircularArray *ca, U64 value) ;
void circularScan(CircularArray *ca, U64 *minimum, size_t *position) ;

// Deque to handle inserted values
Deque *dequeCreate(size_t size) ;
static void dequeDestroy (Deque *dq) { free (dq) ; }
void dequeInsert(Deque *dq, U64 value) ;
bool dequeGetMin(Deque *dq, U64 *minimum, size_t *s_mer_position) ;


Seqhash *seqhashCreate (int k, int w, int seed) ;
static void seqhashDestroy (Seqhash *sh) { free (sh) ; }
  
 // simple iterator to return all kmers
 SeqhashIterator *seqhashIterator (Seqhash *sh, char *s, int len) ;
 bool seqhashNext (SeqhashIterator *si, U64 *kmer, size_t *pos, bool *isF) ;
 
 static void seqhashIteratorDestroy (SeqhashIterator *si)
 { free (si->hash) ; free (si->isForward) ; free (si) ;}
 
 static inline U64 hashRC (SeqhashIterator *si, bool *isForward);
 static inline U64 advanceHashRC (SeqhashIterator *si, bool *isForward);
 // (closed) syncmer extracts w-mers that end with a minimal kmer
 // these provide a cover, and have good distribution properties
 SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len) ;
 bool syncmerNext (SeqhashIterator *si, U64 *kmer, size_t *k_pos, size_t *s_pos, bool *isF) ;

 SeqhashIterator *syncmerIterator_original (Seqhash *sh, char *s, int len) ;
 bool syncmerNext_original (SeqhashIterator *si, U64 *kmer, size_t *pos, bool *isF) ;
 
 SeqhashIterator *syncmerNaiveIterator (Seqhash *sh, char *s, int len) ;
 bool syncmerNaiveNext (SeqhashIterator *si, U64 *smer, size_t *kmer_pos, size_t *smer_pos, bool *isF);

 SeqhashIterator *syncmerDequeIterator (Seqhash *sh, char *s, int len) ;
 bool syncmerDequeNext (SeqhashIterator *si, U64 *smer, size_t *kmer_pos, size_t *smer_pos, bool *isF);

 // utilities
 static inline U64 kHash (Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1) ; }


static inline uint8_t base_to_bits(char base) {
  switch(base) {
      case 'A': case 'a': return 0;
      case 'C': case 'c': return 1;
      case 'G': case 'g': return 2;
      case 'T': case 't': return 3;
      default: return 0; // Treat Ns and  unknown as 'A'
  }
}