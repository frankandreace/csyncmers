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

#include "closed_syncmers_syng.h"

/************** basic hash functions *************/

static inline U64 hashRC (SeqhashIterator *si, bool *isForward)
{ U64 hashF = kHash (si->sh, si->h) ;
  U64 hashR = kHash (si->sh, si->hRC) ;
#ifdef DEBUG
  printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", si->h, si->hRC, hashF, hashR) ;
#endif
  if (hashF < hashR) { *isForward = true ; return hashF ; }
  else { *isForward = false ; return hashR ; }
}

static inline U64 advanceHashRC (SeqhashIterator *si, bool *isForward)
{ Seqhash *sh = si->sh ;
  if (si->s < si->sEnd)
    { si->h = ((si->h << 2) & sh->mask) | base_to_bits(*si->s) ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[(int)base_to_bits(*si->s)] ;
      ++si->s ;
      return hashRC (si, isForward) ;
    }
  else
    return U64MAX ;
}
/************ iterator to run across a sequence, returning (a subset of) hashes ***********/

/*************** this basic one returns all the hashes *********************/
/*************** and its creator is used as a base by the others ***********/

Seqhash *seqhashCreate (int k, int w, int seed)
{
  assert (sizeof (U64) == 8) ;
  Seqhash *sh = new0 (1, Seqhash) ;
  sh->k = k ; if (k < 1 || k >= 32) die ("seqhash k %d must be between 1 and 32\n", k) ;
  sh->w = w ; if (w < 1) die ("seqhash w %d must be positive\n", w) ;
  sh->seed = seed ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  srandom (seed) ;
  sh->factor1 = (random() << 32) | random() | 0x01 ;
  sh->shift1 = 64 - 2*k ;
  sh->factor2 = (random() << 32) | random() | 0x01 ;
  sh->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { sh->patternRC[i] = (3-i) ; sh->patternRC[i] <<= 2*(k-1) ; }
  return sh ;
}

/************ iterators to run across a sequence, returning (a subset of) hashes ***********/

/*************** this basic one returns all the hashes *********************/
/*************** and its creator is used as a base by the others ***********/

SeqhashIterator *seqhashIterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashIterator *si = new0 (1, SeqhashIterator) ;
  si->sh = sh ;
  si->s = s ; si->sEnd = s + len ;
  si->hash = new0 (sh->w, U64) ;
  si->isForward = new0 (sh->w, bool) ;
  int count = 0; // DEBUG
  if (len < sh->k)
    si->isDone = true ; // edge case
  else
    { int i ;			/* preinitialise the hashes for the first kmer */
      for (i = 0 ; i < sh->k ; ++i, ++si->s)
	{ si->h = (si->h << 2) | *si->s ;
	  si->hRC = (si->hRC >> 2) | sh->patternRC[(int)*(si->s)] ; count++;
	}
      *si->hash = hashRC (si, si->isForward) ;
    }
  printf("PRENINITIALIZED FOR %d NUCLEOTIDES.\n", count);
  printf("HASH VALUE IS %llu\n",si->h) ;
  return si ;
}

bool seqhashNext (SeqhashIterator *si, U64 *kmer, size_t *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

  if (kmer) *kmer = *si->hash ; 
  if (pos) *pos = si->iStart ;
  if (isF) *isF = *si->isForward ;

  if (si->s >= si->sEnd)
    si->isDone = true ;
  else
    { *si->hash = advanceHashRC (si, si->isForward) ;
      ++si->iStart ;
      printf("NEW HASH IS %llu\n", *si->hash) ;
    }
  
  return true ;
}

/************ same for closed syncmer ***********/

SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len)
{
  SeqhashIterator *si = seqhashIterator (sh, s, len) ;
  if (len < sh->w + sh->k) si->isDone = true ; // because we are looking for (w+k)-mers not k-mers here
  if (si->isDone) return si ;
    
  /* store first w hashes in hash and set ->min */
  si->min = si->hash[0] ;
  printf("First hash is %llu\n", si->min) ;
  int i ;
  int count = 0;
  for (i = 1 ; i < sh->w ; ++i)
    { si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      printf("Hash of %d is %llu\n", i, si->hash[i]) ;
      count++;
      if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }
  fprintf(stderr,"# loops initialization syncmeriterator: %d\n", count);
  count = 0;
  // si->iStart = 0 ; // from initialisation
  if (si->hash[0] == si->min || si->hash[sh->w-1] == si->min)  {printf("NO WHILE\n") ; return si ; }// we are done
  while (true)
    { U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      count++;
      if (si->s >= si->sEnd) { si->isDone = true ; printf("# iterations first syncmer: %d\n", count); return si ; }
      si->hash[si->iStart++] = x ;
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ; printf("# iterations first syncmer: %d\n", count); return si ; }
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
    {
      printf("# iterations first syncmer: %d\n", count);
	return si ;
    }
    }
  die ("syncmer initialisation failure") ;
}

// bool syncmerNext (SeqhashIterator *si, U64 *kmer, size_t *k_pos, size_t *s_pos, bool *isF)
// {
//   if (si->isDone) return false ; /* we are done */

// #ifdef DEBUG
//   printf ("base %d, iStart %d, min %" PRIx64 "\n", si->base, si->iStart, si->min) ;
//   int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
//   printf ("\n") ;
// #endif

//   if (kmer) *kmer = si->hash[si->iStart] ;
//   if (k_pos) *k_pos = si->base + si->iStart ;
//   if (isF) *isF = si->isForward[si->iStart] ;

//   if (si->hash[si->iStart] == si->min) // need to find new min - could use a heap, but not so bad to search here
//     { int i ;
//       si->min = si->hash[si->iStart] = U64MAX ;
//       for (i = 0 ; i < si->sh->w ; ++i) if (si->hash[i] < si->min) si->min = si->hash[i] ;
//     }
  
//   while (true) // move forwards to the next minimum
//     { 
//       U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
//       printf("NEW HASH IS %llu\n", si->h) ;
//       if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
//       si->hash[si->iStart++] = x ;
//       if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }
//       if (x <= si->min) // min at the end of the w-mer
// 	{ si->min = x ;
// 	  	  printf (" syncmerNext   at_end %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
//       if (s_pos) *s_pos = *k_pos ; // the minimal s-mer is at the beginning of the k-mer
// 	  return true ;
// 	}
//       if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
// 	{
// 	  printf (" syncmerNext at_start %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
//       if (s_pos) *s_pos = *k_pos ; // the minimal s-mer is at the beginning of the k-mer
// 	  return true ;
// 	}
//     }
// }

bool syncmerNext (SeqhashIterator *si, U64 *kmer, size_t *k_pos, size_t *s_pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, min %" PRIx64 "\n", si->base, si->iStart, si->min) ;
  int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
  printf ("\n") ;
#endif

  if (kmer) *kmer = si->hash[si->iStart] ;
  if (k_pos) *k_pos = si->base + si->iStart ;
  if (isF) *isF = si->isForward[si->iStart] ;

  if (si->hash[si->iStart] == si->min) // need to find new min - could use a heap, but not so bad to search here
    { int i ;
      si->min = si->hash[si->iStart] = U64MAX ;
      for (i = 0 ; i < si->sh->w ; ++i) if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }
  
  while (true) // move forwards to the next minimum
    { U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
      si->hash[si->iStart++] = x ;
      if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ;
	  //	  printf (" syncmerNext   at_end %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
	  if (s_pos) *s_pos = *k_pos ; // the minimal s-mer is at the beginning of the k-mer
    return true ;
	}
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	{
	  // printf (" syncmerNext at_start %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
    if (s_pos) *s_pos = *k_pos ; // the minimal s-mer is at the beginning of the k-mer
	  return true ;
	}
    }
}


void compute_closed_syncmers_syng(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results) {
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }
    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    // U64 length = (uint64_t)len;
    *num_results = 0;
    size_t window_size = (U64)K - (U64)S + 1;

    // Durbin's function works with the parameters as follows:
    // its K is normal S
    // its W is the window size, i.e. K-S (probably 0-based w+1 elements)

    //hash structure (which initializes k, w, mask, etc.)
    // fprintf(stderr,"The seed is %llu\n", seed);
    // fprintf(stderr,"The k-mer length is %d\n", K);
    // fprintf(stderr,"The s-mer length is %d\n", S);
    // fprintf(stderr,"The syng window size is %lu\n", window_size);

    Seqhash *sh = seqhashCreate(S, window_size, seed);

    fprintf(stderr,"Seqhash has been created\n");
    // initializing the syncmer iterator
    SeqhashIterator *si = syncmerIterator(sh, sequence_input, len);

    fprintf(stderr,"Syncmer Iterator has been created\n");
    U64 kmer ;
    size_t k_pos ;
    size_t s_pos ;
    bool isF ;
    size_t count = 0 ;

    while (syncmerNext(si, &kmer, &k_pos, &s_pos, &isF))
    {
        count++;
        printf("KMER: %llu, K_POS: %lu\n", kmer, k_pos) ;
        add_minimizer(results, num_results, kmer, k_pos, s_pos);
    }
    seqhashIteratorDestroy(si);
    seqhashDestroy(sh);

    for(int i = 0; i < count; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
   }

}

void compute_closed_syncmers_deque(char *sequence_input, int len, int K, int S, MinimizerResult *results, int *num_results) {
    if(len < K) {
        fprintf(stderr, "Sequence length is less than K\n");
        return;
    }

    // setting the seed to 7 as in Durbin's
    U64 seed  = 7;
    *num_results = 0;

    fprintf(stderr,"RUNNING DEQUE\n");

    size_t num_s_mers = len - S + 1;
    U64 *s_mer_hashes = (U64 *)malloc(num_s_mers * sizeof(U64));
    size_t window_size = K - S + 1;

    Seqhash *sh = seqhashCreate(S, window_size, seed);
    SeqhashIterator *si = seqhashIterator(sh, sequence_input, len);

    U64 smer;
    size_t s_pos;
    bool isF;
    int ii = 0;

    // HASH ALL S-MERS
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
            printf("ITERATION %lu: min_pos: %lu; min_val: %llu; start_pos: %lu, start_val: %llu; end_pos: %lu, end_val: %llu\n", kmer_pos, min_pos, s_mer_hashes[min_pos], kmer_pos, s_mer_hashes[kmer_pos], kmer_pos + K - S, s_mer_hashes[kmer_pos + K - S]) ;
            if (min_pos == kmer_pos || min_pos == kmer_pos + K - S){
                add_minimizer(results, num_results, s_mer_hashes[min_pos], kmer_pos, min_pos);
            }
        }
    }

    for(int i = 0; i < num_s_mers; i++){
        printf("%llu\n", s_mer_hashes[i]);
    }
    for(int i = 0; i < *num_results; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
    }

    free(s_mer_hashes);
    free(deque);
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
    printf("RUNNING NAIVE\n");


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
        s_mer_hashes[s_pos] = smer;
    }
    printf("Tot number of iterations: %lu; num_s_mers: %lu\n", s_pos+1, num_s_mers);

    size_t front = 0, back = 0;
    size_t min_pos;
    U64 min_smer;

    // USE ARRAY SCAN TO COMPUTE SYNCMERS IN O(N*K)
    printf("Num kmers: %lu\n", num_k_mers) ;
    printf("WINDOW SIZE IS %lu\n", window_size) ;
    for(size_t i = 0; i < num_k_mers; i++){
        min_smer = U64MAX;
        for (size_t j = i; j < i + window_size; j++){
            printf("%lu,", j) ;
            if (s_mer_hashes[j] < min_smer){
                min_pos = j;
                min_smer = s_mer_hashes[j];
            }
        }
        printf("\n");
        printf("ITERATION %lu: min_pos: %lu; min_val: %llu; start_pos: %lu, start_val: %llu; end_pos: %lu, end_val: %llu\n", i, min_pos, min_smer, i, s_mer_hashes[i], i+K, s_mer_hashes[i+window_size]) ;
        if (min_pos == i || min_pos == i + window_size-1){
            add_minimizer(results,num_results,min_smer,i,min_pos);
        } 
    }

    // DEBUG PRINTS
    for(int i = 0; i < num_s_mers; i++){
        printf("%llu\n", s_mer_hashes[i]);
    }
    for(int i = 0; i < *num_results; i++){
        printf("%llu\t%lu\t%lu\n", results[i].minimizer_hash, results[i].kmer_position, results[i].smer_position);
    }

    free(s_mer_hashes) ;
}