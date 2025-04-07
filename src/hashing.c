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

/************ circular array to handle operations on iterators ***********/
CircularArray *circularArrayCreate(size_t size){
  CircularArray *ca = new0 (1, CircularArray) ;
  ca->hashVector = new0 (size, U64) ;
  ca->window_size = size ;
  ca->current_position = 0 ;
  return ca ;
}

void circularInsert(CircularArray *ca, U64 value){
  ca->hashVector[ca->current_position] = value;
  // printf("IS VALUE (%llu) in HASVECTOR[%lu] (%llu)?\n",value ,ca->current_position, ca->hashVector[ca->current_position]) ;
  ca->current_position++ ;
  if (ca->current_position == ca->window_size) ca->current_position = 0;
}

void circularScan(CircularArray *ca, U64 *minimum, size_t *position){
  U64 current_minimum = ca->hashVector[ca->current_position] ;
  size_t current_minimum_position = 0 ; 
  U64 scan_position ;
  // printf("[CIRCULARSCAN] POS %lu, VAL %llu\n",current_minimum_position,  current_minimum) ;
  for (int i = 1 ; i < ca->window_size ; i++ ) {
    scan_position = (i + ca->current_position) % ca->window_size ;
    // printf("SCAN POSITION IS %lu", scan_position) ;
    // printf("[CIRCULARSCAN] POS %lu, VAL %llu\n",i,  ca->hashVector[scan_position]) ;
    if ( ca->hashVector[scan_position] < current_minimum ){
      current_minimum = ca->hashVector[scan_position] ;
      current_minimum_position = i;
    }
  }
  if (minimum) *minimum = current_minimum ;
  if (position) *position = current_minimum_position ;
}

/************ deque to handle operations on iterators ***********/
Deque *dequeCreate(size_t size) {
  Deque *dq = new0 (1, Deque) ;
  dq->hashVector = new0 (size, U64) ;
  dq->idVector = new0 (size, size_t) ;
  dq->front = 0 ;
  dq->back = 0 ;
  dq->smer_pos = 0;
  dq->window_size = size ;
  return dq ;
}

void dequeInsert(Deque *dq, U64 value) {
  // remove minimum if out of context
  if (dq->smer_pos >= dq->window_size && dq->idVector[dq->front % dq->window_size] <= dq->smer_pos - dq->window_size)
  {
    // printf("MOVING FRONT++\n") ;
    dq->front++;
  }

  // remove elements greater than value
  // printf("DQ->BACK %lu ; DQ->FRONT %lu ; VALUE: %llu \n", dq->back, dq->front, value) ; //; HV[%lu]: %llu
  // , (dq->back % dq->window_size) -1, dq->hashVector[(dq->back % dq->window_size) -1]
  while(dq->back > dq->front && dq->hashVector[(dq->back - 1) % dq->window_size] > value)
  {
    // printf("DECREASING BACK\n") ;
    dq->back--;
  }

  // insert value
  dq->used_pos = dq->back % dq->window_size ;
  dq->hashVector[dq->used_pos] = value ;
  dq->idVector[dq->used_pos] = dq->smer_pos++ ;
  dq->back++;
  // printf("DEQUE FRONT: %llu, BACK: %llu\n", dq->front%dq->window_size, dq->used_pos);
  // printf("INSERTING VALUE %llu AT %lu ( %lu div %lu)\n", value, dq->used_pos, dq->back, dq->window_size) ;
  // printf("WINDOW SIZE: %lu, SMER_POS: %lu\n", dq->window_size, dq->smer_pos-1);
  // printf("[ ");
  // for(int i = 0; i < dq->back - dq->front; i++){
  //   printf(" (v: %llu , %lu ) ,", dq->hashVector[(dq->front + i)% dq->window_size], dq->idVector[(dq->front +i) % dq->window_size]) ;
  // }
  // printf(" ]\n") ;
}

bool dequeGetMin(Deque *dq, U64 *minimum, size_t *s_mer_position) {
  // printf("DEQUE FRONT IS %lu\n", dq->front%dq->window_size) ;
  // printf("[ ");
  // for(int i = 0; i < dq->back - dq->front; i++){
  //   printf(" (v: %llu , %lu, %lu ) ,", dq->hashVector[(dq->front + i)% dq->window_size], dq->idVector[(dq->front +i) % dq->window_size], (dq->front + i)% dq->window_size) ;
  // }
  // printf(" ]\n") ;
  if (dq->smer_pos - 1 >= dq->window_size -1){
    U64 minimum_position = dq->idVector[dq->front%dq->window_size];
    // printf("DEQUE MIN POS IS %lu, S-MER POS IS %lu, min at queue pos %lu\n", minimum_position, dq->smer_pos, dq->front%dq->window_size) ;
    // U64 kmer_position = position - dq->window_size + 1;
    if (minimum_position == dq->smer_pos - dq->window_size || minimum_position == dq->smer_pos - 1){
      if (minimum)  *minimum = dq->hashVector[dq->front%dq->window_size] ;
      if (s_mer_position) *s_mer_position = minimum_position ;
      return true ;
    }
  }
  return false ;
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
  for (i = 0 ; i < 4 ; ++i) { 
    sh->patternRC[i] = (3-i) ; 
    sh->patternRC[i] <<= 2*(k-1) ; 
  }
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
  si->counter = 0; // DEBUG
  if (len < sh->k)
    si->isDone = true ; // edge case
  else
    { 
      int i ;			/* preinitialise the hashes for the first kmer */
      for (i = 0 ; i < sh->k ; ++i, ++si->s)
	    { 
        
        si->h = (si->h << 2) | base_to_bits(*si->s) ;
	      si->hRC = (si->hRC >> 2) | sh->patternRC[(int)base_to_bits(*si->s)] ;
	    }
      *si->hash = hashRC (si, si->isForward) ;

      si->counter++;
    }

  return si ;
}

bool seqhashNext (SeqhashIterator *si, U64 *kmer, size_t *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

  if (kmer) *kmer = *si->hash ; 
  if (pos) *pos = si->iStart ;
  if (isF) *isF = *si->isForward ;

  if (si->s >= si->sEnd){
    si->isDone = true ;
  }
  else
    { 
      *si->hash = advanceHashRC (si, si->isForward) ;
      si->counter++;
      ++si->iStart ;
    }
  
  return true ;
}

/************ same for closed syncmer ***********/

SeqhashIterator *syncmerIterator (Seqhash *sh, char *s, int len)
{
  SeqhashIterator *si = seqhashIterator (sh, s, len) ;
  if (len < sh->w + sh->k - 1) 
    {
      // printf("LEN IS %d, (w+k) is %lu\n", len, sh->w + sh->k);
      si->isDone = true ; // because we are looking for (w+k)-mers not k-mers here
      si->iMin = U64MAX;
      si->min = U64MAX;
      return si ;
    }
    
  /* store first w hashes in hash and set ->min */
  si->min = si->hash[0] ;
  si->iMin = 0;
  si->abs_k_pos = 0;
  si->base = 0;
  
  // computing the first window
  int i ;
  si->iStart = 1;
  for (i = 1 ; i < sh->w ; ++i)
    { si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      #ifdef DEBUG
      printf("[STARTING_WINDOW] NEW HASH IS %llu, position is: %lu, remaining is:%s\n", si->hash[i], i, si->s) ;
      #endif
      si->iStart++ ;
      if (si->hash[i] < si->min){ 
        // printf("OLD MIN IS: %llu at iMin: %lu; NEW MIN IS: %llu at iMin: %lu\n", si->min, si->iMin, si->hash[i], i) ;
        si->min = si->hash[i] ; 
        si->iMin = i ; 
      }
    }
  // printf("MIN IS: %llu at iMin: %lu\n", si->min, si->iMin) ;
  // printf("CURRENT MIN IS %llu at position %lu\n", si->min, si->iMin) ;

  // I should be at the end of the circular array, reset iStart and base.
  if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }

  // If is syncmer, return
  if (si->iMin == 0 || si->iMin == sh->w - 1)  { return si ; } // we are done

  // No syncmers in the entire sequence (IN CASE IT IS JUST 1 K-MER)
  if (si->s >= si->sEnd) {
    si->iMin = U64MAX;
    si->min = U64MAX;
    si->isDone = true ;  
    return si ; 
  }

  // No syncmer in the 1st window, now check in the next ones until I find one or the sequence ends
  while (true)
    { 
      // compute hash of the next nucleotide
      U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      si->abs_k_pos++;
      si->hash[si->iStart++] = x ;

      #ifdef DEBUG
      printf("[INITIALIZATION, AFTER W] NEW HASH IS %llu, position is: %lu; curr_min: %llu, min_pos: %lu, first:%lu, last:%lu\n", x, sh->w - 1 + si->abs_k_pos, si->min, si->iMin , si->base + si->iStart - si->sh->w, sh->w - 1 + si->abs_k_pos ) ;
      #endif

      // if iStart gets to the end (w), reset it to 0 and add w to base.
      if (si->iStart == si->sh->w) { 
        si->base += si->sh->w ; si->iStart = 0 ; 
      }

      // min at the end of the w-mer: it is a syncmer
      if (x < si->min)
	    { 
        // update min and return
        si->iMin = si->base + si->iStart - 1;
        si->min = x ; 
        // printf("IS SYNCMER. NEW MIN: %llu, position: %lu\n", si->min, si->iMin) ;
        return si ; 
      }

      // min at the beginning of the w-mer: it is a syncmer
      if (si->base + si->iStart - si->sh->w == si->iMin) 
      {
        // return
        // printf("IS SYNCMER. NEW MIN: %llu, position: %lu", si->min, si->iMin) ;
	      return si ;
      }

      // No syncmers in the sequence
      if (si->s >= si->sEnd) {
        si->iMin = U64MAX;
        si->min = U64MAX;
        si->isDone = true ;  
        return si ; 
      }
    }
  die ("syncmer initialisation failure") ;
}

bool syncmerNext (SeqhashIterator *si, U64 *smer, size_t *kmer_pos, size_t *smer_pos, bool *isF)
{
  // update s-mer value, s-mer value position and k-mer position
  if (smer) *smer = si->min ;
  if (kmer_pos) *kmer_pos = si->abs_k_pos;
  if (smer_pos) *smer_pos = si->iMin ;

  // if there is nothing more to read, return
  if (si->s >= si->sEnd) { 
    si->isDone = true ; 
    return false ; 
 }
  // current min is out of scope (out of the window)
  // need to find new min - could use a heap, but not so bad to search here
  if (si->iStart + si->base == si->iMin + si->sh->w)
    { 
      int i ;
      int looking_position ;
      size_t curr_min ;

      si->min = si->hash[si->iStart] = U64MAX ;
      si->iMin = si->iMin + si->sh->w ;

      // scanning the circular array for a new min
      for (i = 1 ; i < si->sh->w; ++i){

        // start from the 1st value after the one just erased
        looking_position = (si->iStart + i) % si->sh->w;

        // update the new min
        if (si->hash[looking_position] < si->min){
          curr_min = si->base - si->sh->w + i + si->iStart;
          si->iMin = curr_min ;
          si->min = si->hash[looking_position] ;
        }
      }
    }

  // move forwards to the next minimum
  while (true)
  { 
    // compute hash and update circular array
    U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
    si->hash[si->iStart++] = x ;
    si->abs_k_pos++; // advancing k-mer position

    #ifdef DEBUG
    printf("NEW HASH IS %llu, position is: %lu\n", x, si->iStart + si->base - 1) ;
    #endif

    // If at the end of the circular array, reset iStart and base.
    if (si->iStart == si->sh->w) { si->base += si->sh->w ; si->iStart = 0 ; }

    // min at the end of the w-mer
    if (x < si->min)
    { 
      // update min and min position
      si->min = x ;
      si->iMin = si->base + si->iStart - 1;

      // if (si->s >= si->sEnd) { si->isDone = true ;}
      return true ;
    }

    // min at the beginning of the w-mer
    if (si->base + si->iStart - si->sh->w == si->iMin)
    {
      // if (si->s >= si->sEnd) { si->isDone = true ;}
      return true ;
    }

    // if at the end of the sequence and the last computed is not a syncmer, 
    // return false to stop after yielding the last correct syncmer pos

    if (si->s >= si->sEnd) { 
      si->isDone = true ;
      return false ; 
   }

  }
}


SeqhashIterator *syncmerNaiveIterator (Seqhash *sh, char *s, int len){
  
  // initialize seqhashiterator
  SeqhashIterator *si = seqhashIterator(sh,s,len) ;
  si->ca = circularArrayCreate(si->sh->w) ;
  circularInsert(si->ca, si->hash[0]) ;
  si->abs_k_pos = 0 ;

  // Compute hashes for the 1st window
  for (int i = 1; i < si->sh->w ; i++ ){
    U64 x = advanceHashRC(si, si->isForward) ;
    // printf("HASH COMPUTED IS %llu\n", x) ;
    circularInsert(si->ca, x) ;
  }

  // verify syncmer condition
  circularScan(si->ca, &si->min, &si->iMin) ;
  // printf("[1st scan] MIN IS %llu, at %lu\n", si->min, si->iMin) ;
  if (si->iMin == 0 || si->iMin == si->sh->w - 1){
    return si;
  }

  // if end of the sequence, return
  if (si->s >= si->sEnd) {
    si->iMin = U64MAX;
    si->min = U64MAX;
    si->isDone = true ;  
    return si ; 
  }

  // if not syncmer, compute next hash until syncmer found or end of sequence
  while(true){
    U64 x = advanceHashRC(si, si->isForward) ;
    // printf("HASH COMPUTED IS %llu\n", x) ;
    circularInsert(si->ca, x) ;
    circularScan(si->ca, &si->min, &si->iMin) ;
    si->abs_k_pos++ ;
    // if syncmer, return
    // printf("[WHILE INITIALIZATION] MIN IS %llu, at %lu\n", si->min, si->iMin) ;
    if (si->iMin == 0 || si->iMin == si->sh->w - 1){
      si->iMin = si->iMin + si->abs_k_pos ;
      // printf("RETURNING SI\n") ;
      return si;
    }
    // if end of sequence, stop
    if (si->s >= si->sEnd) {
      si->iMin = U64MAX;
      si->min = U64MAX;
      si->isDone = true ;  
      return si ; 
    }
  }
}

bool syncmerNaiveNext(SeqhashIterator *si, U64 *smer, size_t *kmer_pos, size_t *smer_pos, bool *isF){
  // update s-mer value, s-mer value position and k-mer position
  if (smer) *smer = si->min ;
  if (kmer_pos) *kmer_pos = si->abs_k_pos;
  if (smer_pos) *smer_pos = si->iMin ;
  // printf("UPDATED. MIN IS %llu, at %lu, kpos: %lu\n", si->min, si->iMin, si->abs_k_pos) ;
  // if there is nothing more to read, return
  if (si->s >= si->sEnd) { 
    si->isDone = true ; 
    return false ; 
  }

  // keep scanning for the next syncmer
  while(true){
    U64 x = advanceHashRC(si, si->isForward) ;
    circularInsert(si->ca, x) ;
    circularScan(si->ca, &si->min, &si->iMin) ;
    si->abs_k_pos++;

    // if syncmer, return
    if (si->iMin == 0 || si->iMin == si->sh->w - 1){
      si->iMin = si->abs_k_pos + si->iMin; // update the syncmer position
      // printf("RETURNIG. FOUND MIN %llu, at %lu, kpos: %lu\n", si->min, si->iMin, si->abs_k_pos) ;
      return true;
    }

    // if end of sequence, and no final syncmer stop
    if (si->s >= si->sEnd) {
      si->iMin = U64MAX;
      si->min = U64MAX;
      si->isDone = true ;  
      return false ; 
    }
  }

}

SeqhashIterator *syncmerDequeIterator (Seqhash *sh, char *s, int len) {
  // initialize seqhashiterator
  // printf("SEQHASHITERATOR\n") ;
  SeqhashIterator *si = seqhashIterator(sh,s,len) ;
  // printf("DEQUE CREATE\n") ;
  si->dq = dequeCreate(si->sh->w) ;
  si->abs_k_pos = 0 ;

  dequeInsert(si->dq, si->hash[0]) ;

  // Compute all hashes except last for the 1st window 
  for (int i = 1; i < si->sh->w - 1; i++ ){
    U64 x = advanceHashRC(si, si->isForward) ;
    // printf("HASH COMPUTED IS %llu\n", x) ;
    dequeInsert(si->dq, x) ;
  }

  // until it is not yielding a syncmer
  // printf("FIRST WHILE\n") ;
  while(true){
    U64 x = advanceHashRC(si, si->isForward) ;
    // printf("HASH COMPUTED IS %llu\n", x) ;
    dequeInsert(si->dq, x);
    if (dequeGetMin(si->dq, &si->min, &si->iMin)){
        // printf("\nRETURNING MIN %llu at pos %lu\n\n", si->min, si->iMin) ;
        break ;
    }
    si->abs_k_pos++;
    if (si->s >= si->sEnd) {
      // printf("END OF SEQUENCE IN ITIALIZATION\n") ;
      si->iMin = U64MAX;
      si->min = U64MAX;
      si->isDone = true ;  
      return si ; 
    }
  }
  return si ;
}

bool syncmerDequeNext (SeqhashIterator *si, U64 *smer, size_t *kmer_pos, size_t *smer_pos, bool *isF) {
  // update s-mer value, s-mer value position and k-mer position
  if (smer) *smer = si->min ;
  if (kmer_pos) *kmer_pos = si->abs_k_pos;
  if (smer_pos) *smer_pos = si->iMin ;
  // printf("UPDATED. MIN IS %llu, at %lu, kpos: %lu\n", si->min, si->iMin, si->abs_k_pos) ;
  // if there is nothing more to read, return
  if (si->s >= si->sEnd) { 
    si->isDone = true ; 
    return false ; 
  }

  // keep scanning for the next syncmer
  while(true){
    U64 x = advanceHashRC(si, si->isForward) ;
    // printf("HASH COMPUTED IS %llu\n", x) ;
    dequeInsert(si->dq, x);
    si->abs_k_pos++;
    if (dequeGetMin(si->dq, &si->min, &si->iMin)){
      // printf("\nRETURNING MIN %llu at k_pos %lu, s_pos %lu\n\n", si->min, si->abs_k_pos, si->iMin) ;
      return true;
    }
    if (si->s >= si->sEnd) {
      si->iMin = U64MAX;
      si->min = U64MAX;
      si->isDone = true ;  
      return false ; 
    }
  }
}

SeqhashIterator *syncmerIterator_original (Seqhash *sh, char *s, int len)
{
  SeqhashIterator *si = seqhashIterator (sh, s, len) ;
  if (len < sh->w + sh->k) si->isDone = true ; // because we are looking for w-mers not k-mers here
  if (si->isDone) return si ;
    
  /* store first w hashes in hash and set ->min */
  si->min = si->hash[0] ;
  int i ;
  for (i = 1 ; i < sh->w ; ++i)
    { si->hash[i] = advanceHashRC (si, &si->isForward[i]) ;
      if (si->hash[i] < si->min) si->min = si->hash[i] ;
    }

  // si->iStart = 0 ; // from initialisation
  if (si->hash[0] == si->min || si->hash[sh->w-1] == si->min) return si ; // we are done
  while (true)
    { U64 x = advanceHashRC (si, &si->isForward[si->iStart]) ;
      if (si->s >= si->sEnd) { si->isDone = true ; return si ; }
      si->hash[si->iStart++] = x ;
      if (x <= si->min) // min at the end of the w-mer
	{ si->min = x ; return si ; }
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	return si ;
    }
  die ("syncmer initialisation failure") ;
}

bool syncmerNext_original (SeqhashIterator *si, U64 *kmer, size_t *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, min %" PRIx64 "\n", si->base, si->iStart, si->min) ;
  int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hash[j]) ;
  printf ("\n") ;
#endif

  if (kmer) *kmer = si->hash[si->iStart] ;
  if (pos) *pos = si->base + si->iStart ;
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
	  return true ;
	}
      if (si->hash[si->iStart] == si->min) // min at the beginning of the w-mer
	{
	  // printf (" syncmerNext at_start %" PRId64 " %" PRIx64 " %" PRIu64 "\n", si->iStart, x, si->sEnd-si->s) ;
	  return true ;
	}
    }
}