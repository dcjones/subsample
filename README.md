
A small, Unix-y, program to efficiency subsample lines from a file
without replacement.

# Synopsis

`subsample -n 1000 in_file [in_file2...] > out_file`

This produces an output file with 1000 lines sampled uniformly without
replacement from the concatenated input files. Though the selected lines are
random, their order is preserved.

The goal of this program is to draw large samples from huge input files. Unlike
reservoir sampling, or more naive methods, it avoids keeping any significant
portion of the input in memory, requiring just m/8 bytes.  E.g. to sample from a
file with 100 million lines, you need less than 13MB, regardless of the size in
bytes of the input file.


