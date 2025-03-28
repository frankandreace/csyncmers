#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <stdio.h>
#include <stdlib.h>

#define STREAM_BUFSIZ (1<<12)

/* Stream type definitions */
typedef enum {
    STREAM_FILE,
    STREAM_FASTA
} stream_type;

typedef struct stream {
    stream_type type;
    int reverse;      // Not used in this simple version; reserved for future use
    long len;         // Total length of the file
    long pos;         // Current position in the file
    FILE *f;          // File pointer
    long buf_ptr;     // File offset corresponding to the start of the buffer
    size_t buf_size;  // Number of valid bytes in the buffer
    unsigned char buf[STREAM_BUFSIZ];  // Internal buffer
} stream;

/* Opens a FASTA file stream. Returns NULL on failure. */
stream *stream_open_fasta(FILE *fin);

/* Returns the next nucleotide from the FASTA stream:
   - Returns lowercase nucleotide letters ('a', 't', 'g', 'c', 'n').
   - Returns 0 when a new FASTA record is reached (i.e. end of a sequence).
   - Returns EOF when the file is done.
*/
int stream_getnext(stream *S);

/* reads an entire sequence from the fasta file*/
char *read_sequence(struct stream *S);

/* Closes the FASTA stream and cleans up resources. */
void stream_close(stream *S);

#endif /* FASTA_READER_H */
