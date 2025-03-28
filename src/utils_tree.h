#ifndef UTILS_H
#define UTILS_H

#ifndef NEW
 #define NEW(p,n) {p = malloc((n)*sizeof(*p));if ((p)==NULL) {printf("not enough memory\n"); exit(1);};}
#endif

typedef unsigned char uchar;
typedef         short  i16;
typedef unsigned short u16;
typedef          int   i32;
typedef unsigned int   u32;
typedef          long  i64;
typedef unsigned long  u64;

#define Cdollar 0
#define CA 1
#define CC 2
#define CG 3
#define CT 4
#define CN 5
//#define Cgap 6
#ifndef ALPHABET_SIZE
 #define ALPHABET_SIZE 5
#endif
#define Cgap (ALPHABET_SIZE+1)

/////////////////////////////////////////////////////
// $: 文字列の一番最後にだけ現れる (0)
// a, c, g, t, n: 通常文字 (1,..., ALPHABET_SIZE)
// -: 各文字列の先頭を指す (ALPHABET_SIZE+1)
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////
// binary encoding の場合
// $: 文字列の一番最後にだけ現れる (0)
// 0, 1, 2: 通常文字 (1,..., ALPHABET_SIZE=3)
// -: 各文字列の先頭を指す (ALPHABET_SIZE+1)
/////////////////////////////////////////////////////


static int convertchar(uchar t)
{
  int c = 0;
  switch (t) {
  case '$':  c = Cdollar;  break;
  case 'a':
  case 'A':  c = CA;  break;
  case 'c':
  case 'C':  c = CC;  break;
  case 'g':
  case 'G':  c = CG;  break;
  case 't':
  case 'T':  c = CT;  break;
  case 'n':
  case 'N':  c = CN;  break;
  case '-':  c = Cgap;  break;
  default:  //printf("error char = %c [%02x]\n",t,t);
             c = -1;  break;
  }
  return c;
}

#define STREAM_MEM 0
#define STREAM_MMAP 1
#define STREAM_FILE 2
#define STREAM_FASTA 3
typedef struct stream {
  int type;
  int reverse;
  long len;
  long pos;
  union {
    struct {
      uchar *p;
    } mem;
    struct {
      uchar *p;
    } mmap;
    struct {
      FILE *f;
      long ptr;
      uchar *buf;
    } file;
  } u;
}* stream;

static stream stream_open_mem(uchar *p, int reverse)
{
  stream S;
  NEW(S, 1);
  S->type = STREAM_MEM;
  S->reverse = reverse;
  S->u.mem.p = p;
  long len = 0;
  while (*p != 0) {
    p++;
    len++;
  }
  S->len = len;
  if (reverse) {
    S->pos = len-1;
  } else {
    S->pos = 0;
  }
  return S;
}

static int stream_getnext_mem(stream S)
{
  int c = S->u.mem.p[S->pos];
  if (S->reverse) S->pos--; else S->pos++;
  return c;
}


#define STREAM_BUFSIZ (1<<12)
static stream stream_open_file(FILE *fin, int reverse)
{
  stream S;
  NEW(S, 1);
  S->type = STREAM_FILE;
  S->reverse = reverse;
  fseek(fin, 0, SEEK_END);
  long len = ftell(fin);
  fseek(fin, 0, SEEK_SET);
  NEW(S->u.file.buf, STREAM_BUFSIZ);

  S->len = len;
  if (reverse) {
    S->pos = len-1;
  } else {
    S->pos = 0;
  }
  S->u.file.ptr = -1;
  S->u.file.f = fin;
  return S;
}


static int stream_getnext_file(stream S)
{
  if (S->reverse == 0) {
    if (S->u.file.ptr == -1 || S->pos >= S->u.file.ptr + STREAM_BUFSIZ) {
      long ptr = (S->pos / STREAM_BUFSIZ) * STREAM_BUFSIZ;
      fseek(S->u.file.f, ptr, SEEK_SET);
      size_t size = fread(S->u.file.buf, 1, STREAM_BUFSIZ, S->u.file.f);
      if (size == 0) {
        S->u.file.ptr = -1;
        return EOF;
      }
      S->u.file.ptr = ptr;
    }
  } else {
    if (S->u.file.ptr == -1 || S->pos < S->u.file.ptr) {
      long ptr = (S->pos / STREAM_BUFSIZ) * STREAM_BUFSIZ;
      fseek(S->u.file.f, ptr, SEEK_SET);
      size_t size = fread(S->u.file.buf, 1, STREAM_BUFSIZ, S->u.file.f);
      if (size == 0) {
        S->u.file.ptr = -1;
        return EOF;
      }
      S->u.file.ptr = ptr;
    }
  }
  int c = S->u.file.buf[S->pos - S->u.file.ptr];
  if (S->reverse) S->pos--; else S->pos++;
  return c;
}

static void stream_ungetc_file(stream S)
{
  if (S->reverse) S->pos++; else S->pos--;
}

static int stream_getnext_fasta_reverse(stream S)
{
  return EOF;  
}

static int stream_getnext_fasta(stream S)
{
  if (S->reverse) return stream_getnext_fasta_reverse(S);
  int c = EOF;
  while (1) {
    c = stream_getnext_file(S);
    if (c == EOF) break;
    if (c == '\n') {
      c = stream_getnext_file(S);
      if (c == EOF) return c;
      if (c == '>') {
        stream_ungetc_file(S);
        return 0; // 1つの文字列の最後なら 0 を返す．
      }
    }
    if (c == '>') { // コメントの先頭
      while (1) { // 改行まで読み飛ばす
        c = stream_getnext_file(S);
        if (c == EOF) return c;
        if (c == '\n') break;
      }
    }
    c = tolower(c);
    if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n') break;
  }
  return c;
}



static int stream_getnext(stream S)
{
  if (S->pos < 0 || S->pos >= S->len) return EOF;
  switch (S->type) {
    case STREAM_MEM: return stream_getnext_mem(S);
    case STREAM_FILE: return stream_getnext_file(S);
    case STREAM_FASTA: return stream_getnext_fasta(S);
    default: printf("stream_getnext: unknown type %d\n", S->type);
             exit(1);
  }
  return EOF;
}

static void stream_close(stream S)
{
  switch (S->type) {
    case STREAM_MEM:
      break;
    case STREAM_FILE:
    case STREAM_FASTA:
      free(S->u.file.buf);
      //fclose(S->u.file.f);
      break;
    default: printf("stream_close: unknown type %d\n", S->type);
             exit(1);
  }
  free(S);
}

#undef STREAM_BUFSIZ

#endif
