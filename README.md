# CSyncmers

Comparison of computation speed of closed syncmers in pure C. Supports arbitrary k values and s < 64 (using Syng hash function). Implementations are integrated into syng's hashing fucntions and structures.

Rescan is faster than deque.

#### Folder
- Executables are in bin/;
- Test to check fucntions is test.sh;
- Source files are in src/;
- Data is in data/;


#### How to run
To test correctness do:

```
Make all
```
This creates (in bin/) a test executable that tests the correctness of syncmer computation for given K, S and sequence.
In test.sh you can change the K and S value, as well as the randoms sequence length.

200 test of various K and S are performed.

If you do:

```
./test.sh 
```

To test speed do:

```
Make clean
Make sync
```
This creates (in bin/) a sync executable that test speed of syncmer computation for given K, S and fasta file.
In test_speed you can change the file you want to test on.

and then do

```
./speed_test.sh 
```

### Some of existing code

Strobealign code: uses 64-bit smers 
https://github.com/ksahlin/strobealign/blob/71866c31b578e5166c83aaf1fde79d238246490d/src/randstrobes.cpp#L57

Minimap2 code: uses 64-bit smers
https://github.com/lh3/minimap2/blob/c2f07ff2ac8bdc5c6768e63191e614ea9012bd5d/sketch.c#L145-L192

Curiouscoding blog:
https://curiouscoding.nl/posts/fast-minimizers/

Sliding window minimum algorithm explanation:
https://github.com/keegancsmith/Sliding-Window-Minimum?tab=readme-ov-file#sliding-window-minimum-algorithm

SYNG: 
https://github.com/richarddurbin/syng/

### TODOS:
- [x] Integrate Syng rolling hash function into Rayan's code
- [x] Test Syng syncmer extraction
- [x] Fix syncmer iterator as it returns 0 reverse complement hash values.
- [x] Read curiouscoding blog and watch curiouscoding video
- [x] Integrate deque, rescan and naive to Syng hashing and test speed
