#include <time.h>
#include <stdio.h>

void print_timing_stats(char* name, double elapsed_time){
    printf("Processing time: %.2f seconds\n", elapsed_time);
}


void print_troughput(){
    //   double file_size_mb = file_size / (1024.0 * 1024.0);
    //   double processing_speed = file_size_mb / elapsed_time;
    //   printf("File size: %.2f MB\n", file_size_mb);
    //   printf("Processing speed: %.2f MB/sec\n", processing_speed);
}

void print_benchmark(char* name, double start_time, double end_time, int syncmer_count, FILE *file){
    
    // precomputing useful variables
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // file size

    printf("%s Performance metrics.\n", name);

    printf("Syncmers found: %d\n", syncmer_count);
    // print time stats
    print_timing_stats(name, elapsed_time);

    // estimated throughput

    printf("\n");
}

//   
// Get file size in bytes
// off_t get_file_size(const char *filename) {
//     struct stat st;
//     if (stat(filename, &st) == 0) {
//       return st.st_size;
//     }
//     return -1;
//   }