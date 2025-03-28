#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

#include "utils_tree.h"

#define WINDOW_LEN 55



typedef unsigned long smer_t;

u64 smer_hash(smer_t sm)
{
  sm = sm & 0xffffffffffffffffUL;
  return sm;
}

smer_t get_next_smer(int smer_len, smer_t sm, uchar c)
{
  c = convertchar(c)-1;
  sm = sm << 2;
  sm += c;
  sm = sm & ((1UL<< 2*smer_len)-1);
  return sm;
}

int get_window(int window_len, char buf[], stream s_in)
{
  int end = 0;
  for (int i=0; i<window_len; i++) {
    int c = stream_getnext(s_in);
    if (c == EOF) {
      end = 1;
      c = 'a';
    }
    buf[i] = c;
  }
  return end;
}

int get_next_window(int window_len, char buf[], stream s_in)
{
  int end = 0;
  for (int i=0; i<window_len-1; i++) {
    buf[i] = buf[i+1];
  }
  int c = stream_getnext(s_in);
  if (c == EOF) {
    end = 1;
    c = 'a';
  }
  buf[window_len-1] = c;
  return end;
}


int check_open_syncmer_naive(int smer_len, int window_len, char buf[])
{
  smer_t sm = 0;
  for (int i=0; i<smer_len; i++) {
    sm = get_next_smer(smer_len, sm, buf[i]);
  }
  u64 sm_hash = smer_hash(sm);
  for (int i=smer_len; i<window_len; i++) {
    sm = get_next_smer(smer_len, sm, buf[i]);
    u64 sm_hash2 = smer_hash(sm);
    //printf("i=%d h1 %ld h2 %ld\n", i, sm_hash, sm_hash2);
    if (sm_hash2 < sm_hash) return 0;
  }
  return 1;
}

#if 0
typedef struct {
  int window_len, smer_len;
  char buf[WINDOW_LEN*2];
  u64 hvalue[WINDOW_LEN];
  int pos[WINDOW_LEN];
}* sorted_smers;

sorted_smers get_sorted_smers(int smer_len, int window_len, char buf[])
{
  sorted_smers ss;
  NEW(ss, 1);
  ss->smer_len = smer_len;
  ss->window_len = window_len;

  smer_t sm = 0;
  for (int i=0; i<smer_len; i++) {
    ss->buf[i] = buf[i];
    sm = get_next_smer(smer_len, sm, buf[i]);
  }
  for (int i=0; i<window_len; i++) {
    ss->buf[smer_len + i] = buf[smer_len + i];
    ss->pos[i] = i;
    ss->hvalue[i] = smer_hash(sm);
    sm = get_next_smer(smer_len, sm, buf[smer_len + i]);
  }
  for (int i=window_len; i<window_len*2-smer_len; i++) {
    ss->buf[smer_len + i] = buf[smer_len + i];
  }
  for (int i=0; i<window_len-1; i++) {
    for (int j=i; j<window_len; j++) {
      if (ss->hvalue[j] < ss->hvalue[i]) {
        u64 h = ss->hvalue[j];
        ss->hvalue[j] = ss->hvalue[i];
        ss->hvalue[i] = h;
        int t = ss->pos[j];
        ss->pos[j] = ss->pos[i];
        ss->pos[i] = t;
      }
    }
  }
  return ss;
}

int find_open_syncmer(int smer_len, int window_len, sorted_smers ss0, sorted_smers ss1, int ofs)
{
  int j0 = 0, j1 = 0;
  //printf("ss0 %s\n", ss0->buf);
  for (int i=0; i<window_len; i++) {
    //printf("i=%d pos = %d hv = %ld\n", i, ss0->pos[i], ss0->hvalue[i]);
  }
  //printf("ss1 %s\n", ss1->buf);
  for (int i=0; i<window_len; i++) {
    //printf("i=%d pos = %d hv = %ld\n", i, ss1->pos[i], ss1->hvalue[i]);
  }


  for (int i=0; i<window_len; i++) {
    u64 h0, h1;
    int p0, p1;
    while (j0 < window_len && ss0->pos[j0] < i) { // window [i,i+w-1]
      j0++;
    }
    p0 = ss0->pos[j0];
    h0 = ss0->hvalue[j0];
    //printf("i=%d j0=%d h0=%ld p0=%d\n", i, j0, h0, p0);
    if (i > 0) {
      for (j1 = 0; j1 < window_len; j1++) {
        //printf("j1=%d h=%ld p=%d\n", j1, ss1->hvalue[j1], ss1->pos[j1]);
        if (ss1->pos[j1] <= i - smer_len) {
          h1 = ss1->hvalue[j1];
          p1 = ss1->pos[j1]+window_len;
          //printf("h1=%ld p1=%d\n", h1, p1);
          break;
        }
        if (ss1->hvalue[j1] >= h0) {
          j1 = window_len;
          break;
        }
      }
      if (j1 < window_len) {
        if (h1 < h0) {
          p0 = p1;
          h0 = h1;
        }
      }
    }
    if (p0 == i) {
      printf("i=%d syncmer (%d) ", i+ofs, smer_len);
      for (int j=0; j<window_len; j++) printf("%c", ss0->buf[i+j]);
      printf("\n");
    }
  }
}
#endif

typedef struct {
  int window_len, smer_len;
  int stack_len, p;
  char buf[WINDOW_LEN*2];
  u64 hvalue[WINDOW_LEN];
  int pos[WINDOW_LEN];
}* smer_list;

smer_list smer_list_new(int window_len, int smer_len)
{
  smer_list ss;
  NEW(ss, 1);
  ss->window_len = window_len;
  ss->smer_len = smer_len;
  ss->stack_len = 0;
  ss->p = 0;
  return ss;
}

void smer_list_push(smer_list ss, int pos, u64 hvalue)
{
  int i = ss->stack_len-1;
  while (i >= 0) {
    if (hvalue >= ss->hvalue[i]) break;
    i--;
  }
  ss->hvalue[i+1] = hvalue;
  ss->pos[i+1] = pos;
  ss->stack_len = i+2;
  ss->p++;
}


smer_list create_smer_list(int smer_len, int window_len, char buf[])
{
  smer_list ss = smer_list_new(window_len,smer_len);

  smer_t sm = 0;
  for (int i=0; i<smer_len; i++) {
    sm = get_next_smer(smer_len, sm, buf[i]);
    ss->buf[i] = buf[i];
  }

  for (int i=0; i<window_len; i++) {
    smer_list_push(ss, i, sm);
    sm = get_next_smer(smer_len, sm, buf[smer_len + i]);
    ss->buf[smer_len + i] = buf[smer_len + i];
    //printf("i=%d ", i);
    //for (int j=0; j<ss->stack_len; j++) {
    //  printf("(%d,%ld)", ss->pos[j], ss->hvalue[j]);
    //}
    //printf("\n");
  }
  for (int i=0; i<window_len-smer_len; i++) {
    ss->buf[smer_len + window_len + i] = buf[smer_len + window_len + i];
  }
  return ss;
}

int find_open_syncmer2(int smer_len, int window_len, char buf[], smer_list ss, int ofs)
{
  int syncmer_count = 0;
  smer_t sm = 0;
  for (int i=0; i<smer_len-1; i++) {
    sm = get_next_smer(smer_len, sm, buf[i]); // hash value in the right window
  }
  smer_t min_sm = 1UL<< 2*smer_len;
  int min_pos = 0;

  int left_min_index = 0; // position in the stack of the left window
  int left_min_pos = ss->pos[left_min_index];
  u64 left_min_hvalue = ss->hvalue[left_min_index];
  int right_min_pos = window_len;
  u64 right_min_hvalue = 1UL<< 2*smer_len;

  for (int i=0; i<window_len; i++) {
    int p = left_min_pos;
    if (right_min_hvalue < left_min_hvalue) p = right_min_pos;
    if (p == i) {
      //printf("i=%d syncmer (%d) ", i+ofs, smer_len);
      //for (int j=0; j<window_len; j++) printf("%c", ss->buf[i+j]);
      //printf("\n");
      syncmer_count++;
    }
    // update the minimum value in the left window
    //if (p <= i) { 
    if (left_min_pos <= i) { 
      left_min_index++;
      if (left_min_index < ss->stack_len) {
        left_min_pos = ss->pos[left_min_index];
        left_min_hvalue = ss->hvalue[left_min_index];  
      }
    }
    // update the minimum value in the right window
    //if (1) {
    //if (i < window_len-smer_len-i+1) {
    if (i < window_len-i) {
        sm = get_next_smer(smer_len, sm, buf[smer_len-1 + i]);
      if (sm < right_min_hvalue) {
        right_min_hvalue = sm;
        right_min_pos = window_len+i;
      }
    } else {
      //right_min_hvalue = 1UL<< 2*smer_len;
    }
    //printf("i=%d left_min_hvalue = %ld left_min_pos = %d\n", i+ofs, left_min_hvalue, left_min_pos);
    //printf("i=%d right_min_hvalue = %ld right_min_pos = %d\n", i+ofs, right_min_hvalue, right_min_pos);
  }
  return syncmer_count;
}


// Get file size in bytes
off_t get_file_size(const char *filename) {
  struct stat st;
  if (stat(filename, &st) == 0) {
    return st.st_size;
  }
  return -1;
}

int main(int argc, char *argv[])
{
  if (argc < 2) {
    printf("Usage: %s <input_fasta_file>\n", argv[0]);
    return 1;
  }

  int window_len = WINDOW_LEN;
  int smer_len = 8;
  char C[] = {'$', 'a', 'c', 'g', 't', 'N', '-'};

  // Get file size for performance measurement
  off_t file_size = get_file_size(argv[1]);
  if (file_size < 0) {
    printf("Error: Cannot determine file size of %s\n", argv[1]);
    return 1;
  }
  
  // Start timing
  clock_t start_time = clock();
  
  FILE *fin = fopen(argv[1], "rb");
  if (!fin) {
    printf("Error: Cannot open file %s\n", argv[1]);
    return 1;
  }
  
  stream s_in = stream_open_fasta(fin, 0);

#if 0
  char buf0[WINDOW_LEN*2];
  get_window(window_len*2, buf0, s_in);
  int i = 0;
  while (1) {
    int f = check_open_syncmer_naive(smer_len, window_len, buf0);
    if (f) {
      printf("i=%d syncmer (%d) %s\n", i, smer_len, buf0);
    }
    int end = get_next_window(window_len*2, buf0, s_in);
    if (end) break;
    i += 1;
    if (i > 50) break;
  }
#endif

#if 0
  char buf[WINDOW_LEN*3];
  sorted_smers ss0=NULL, ss1=NULL;
  get_window(window_len*2, buf, s_in);
  ss1 = get_sorted_smers(smer_len, window_len, buf);
  int ofs = 0;
  while (1) {
    for (int i=0; i<window_len; i++) {
      buf[i] = buf[window_len + i];
    }
    if (ss0 != NULL) free(ss0);
    ss0 = ss1;
    int end = get_window(window_len, buf + window_len, s_in);
    ss1 = get_sorted_smers(smer_len, window_len, buf);

    find_open_syncmer(smer_len, window_len, ss0, ss1, ofs);
    ofs += window_len;
    if (end) break;
    if (ofs > 50) break;
  }
#endif

#if 1
  char buf[WINDOW_LEN*2];
  for (int i=0; i<window_len*2; i++) buf[i] = 0;
  get_window(window_len*2, buf, s_in);
  int ofs = 0;
  int end = 0;
  int syncmer_count = 0;
  
  while (1) {
    smer_list ss = create_smer_list(smer_len, window_len, buf);
    find_open_syncmer2(smer_len, window_len, buf+window_len, ss, ofs);
    if (end) break;
    for (int i=0; i<window_len; i++) {
      buf[i] = buf[window_len + i];
    }
    end = get_window(window_len, buf+window_len, s_in);
    ofs += window_len;
    //if (ofs > 50) break;
    free(ss);
  }
#endif

  stream_close(s_in);
  fclose(fin);
  
  // End timing
  clock_t end_time = clock();
  double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
  
  // Calculate and print performance metrics
  double file_size_mb = file_size / (1024.0 * 1024.0);
  double processing_speed = file_size_mb / elapsed_time;
  
  printf("Performance metrics:\n");
  printf("File size: %.2f MB\n", file_size_mb);
  printf("Processing time: %.2f seconds\n", elapsed_time);
  printf("Processing speed: %.2f MB/sec\n", processing_speed);
  printf("Syncmers found: %d\n", syncmer_count);
  
  // Calculate syncmer density (syncmers per MB)
  double syncmer_density = syncmer_count / file_size_mb;
  printf("Syncmer density: %.2f syncmers/MB\n", syncmer_density);
  
  return 0;
}
