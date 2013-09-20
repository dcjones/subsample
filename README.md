
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

# Usage

`subsample [OPTIONS] [FILE...]`

Sample lines from the concatenation of one or more regular files and print to
standard out.

## Options

`-n N`: Number of chunks to sample. (default: 1)

`-s SEED, --seed=SEED`: Seed the random number generator. By default a fixed
seed is used. Specifying `-seed` with no argument will choose a somewhat unique
seed from the current time and the process id.

`-k K, --chunksize=K`: Sample over groups of this many lines. E.g. if `-k 4` is
used, every four lines is considered an entry. (default: 1)

`-d C, --delimiter=C`: Lines are separated by this character. (default: \n)

`--help`: Print this message.


