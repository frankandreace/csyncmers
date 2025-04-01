#include "hashing.h"

/************** basic hash functions *************/

static inline U64 hashRC (SeqhashIterator *si, bool *isForward)
{ U64 hashF = kHash (si->sh, si->h) ;
  U64 hashR = kHash (si->sh, si->hRC) ;
// #ifdef DEBUG
//   printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", si->h, si->hRC, hashF, hashR) ;
// #endif
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
  printf("HASH VALUE IS %llu\n", *si->hash) ;
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
  // si->iMin = 0;
  // printf("First hash is %llu\n", si->min) ;
  int i ;
  int count = 0;
  for (i = 1 ; i < sh->w ; ++i)
    { si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      printf("Hash of %d is %llu\n", i, si->hash[i]) ;
      count++;
      if (si->hash[i] < si->min){ 
        si->min = si->hash[i] ; 
        // si->iMin == i ; 
      }
    }
  fprintf(stderr,"# loops initialization syncmeriterator: %d\n", count);
  count = 0;
  // si->iStart = 0 ; // from initialisation
  if (si->hash[0] == si->min || si->hash[sh->w-1] == si->min)  {return si ; }// we are done
  while (true)
    { U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      count++;
      if (si->s >= si->sEnd) { si->isDone = true ;  return si ; }
      si->hash[si->iStart++] = x ;
      if (x <= si->min) // min at the end of the w-mer
	    { 
        si->min = x ; return si ; 
      }
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
      {
	      return si ;
      }
    }
  die ("syncmer initialisation failure") ;
}

bool syncmerNext (SeqhashIterator *si, U64 *kmer, size_t *k_pos, bool *isF)
{
  if (si->isDone) return false ; /* End of the string: we are done */

// #ifdef DEBUG
//   printf ("base %d, iStart %d, min %" PRIx64 "\n", si->base, si->iStart, si->min) ;
//   int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
//   printf ("\n") ;
// #endif

  if (kmer) *kmer = si->hash[si->iStart] ;
  if (k_pos) *k_pos = si->base + si->iStart ;
  if (isF) *isF = si->isForward[si->iStart] ;

  if (si->hash[si->iStart] == si->min) // need to find new min - could use a heap, but not so bad to search here
    { int i ;
      si->min = si->hash[si->iStart] = U64MAX ;
      for (i = 0 ; i < si->sh->w ; ++i) if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }
  
  while (true) // move forwards to the next minimum
    { 
      U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      printf("NEW HASH IS %llu\n", x) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
      si->hash[si->iStart++] = x ;
      if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ;
	  	  // printf (" syncmerNext at_end %lu %llu %li\n", si->iStart, x, si->sEnd-si->s) ;
      printf("S_MER at END, last s-mer computed was: %llu, the one at start is %llu\n", x,si->hash[si->iStart]) ;
      printf("s-mer is %llu, the position is %lu.\n", si->min, si->base + si->iStart - 1) ;
      // if (s_pos) *s_pos = *k_pos; // the minimal s-mer is at the end of the k-mer
	  return true ;
	}
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	{
	  // printf (" syncmerNext at_start %lu %llu %li\n", si->iStart, x, si->sEnd-si->s) ;
      printf("S_MER at START, last s-mer computed was: %llu, the one at start is %llu\n", x,si->hash[si->iStart]) ;
      printf("s-mer is %llu, the position is %lu.\n", si->min, si->base + si->iStart - si->sh->w) ;
      // if (s_pos) *s_pos = *k_pos + si->sh->w -1 ; // the minimal s-mer is at the beginning of the k-mer
	  return true ;
	}
    }
}
