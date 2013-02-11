/*
   The MIT License

   Copyright (c) 2013 Daniel C. Jones <dcjones@cs.washington.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/


#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


// Same useful things
// ------------------

void* malloc_or_die(size_t n)
{
    void* p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


void* realloc_or_die(void* p, size_t n)
{
    p = realloc(p, n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


// Random number generation
// ------------------------
// This is the complementary multiply with carry (CMWC) algorithm.

typedef struct rng_t
{
    uint32_t Q[4096];
    uint32_t c;
    uint32_t i;
} rng_t;


void rng_init(rng_t* rng, uint32_t seed)
{
    const uint32_t phi = 0x9e3779b9;
    rng->c = 362436;
    rng->i = 4095;

    rng->Q[0] = seed;
    rng->Q[1] = seed + phi;
    rng->Q[2] = seed + 2 * phi;

    uint32_t i;
    for (i = 3; i < 4096; ++i) {
        rng->Q[i] = rng->Q[i - 3] ^ rng->Q[i - 2] ^ phi ^ i;
    }
}


uint32_t rng_get(rng_t* rng)
{
    uint64_t t, a = UINT64_C(18782);
    uint32_t x, r = 0xfffffffe;
    rng->i = (rng->i + 1) & 4095;
    t = a * rng->Q[rng->i] + rng->c;
    rng->c = t >> 32;
    x = t + rng->c;
    if (x < rng->c) {
        ++x;
        ++rng->c;
    }

    return (rng->Q[rng->i] = r - x);
}


uint32_t rng_uniform_int(rng_t* rng, uint32_t k)
{
    uint32_t scale = UINT32_C(0xffffffff) / k;
    uint32_t r;
    do {
        r = rng_get(rng) / scale;
    } while (r >= k);
    return r;
}


// Random double in [0, 1]
double rng_uniform(rng_t* rng)
{
    return (double) rng_get(rng) / (double) UINT32_MAX;
}


// Random double in [0, 1)
double rng_uniform_pos(rng_t* rng)
{
    return (double) rng_get(rng) / ((double) UINT32_MAX + 1.0);
}


// Random gaussian.
double rng_gaussian(rng_t* rng, double sd)
{
    double x, y, r2;
    do {
        x = -1 + 2 * rng_uniform_pos(rng);
        y = -1 + 2 * rng_uniform_pos(rng);
        r2 = x * x + y * y;
    } while (r2 > 1.0 || r2 == 0.0);
    return sd * y + sqrt(-2.0 * log(r2) / r2);
}


// Random gamma.
double rng_gamma(rng_t* rng, double a, double b)
{
    if (a < 1) {
        double u = rng_uniform_pos(rng);
        return rng_gamma(rng, 1.0 + a, b) * pow(u, 1.0 / a);
    }

    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt(d);

    do {
        do {
            x = rng_gaussian(rng, 1.0);
            v = 1.0 + c * x;
        } while (v <= 0);

        v = v * v * v;
        u = rng_uniform_pos(rng);
    } while (u >= 1 - 0.0331 * x * x * x * x &&
             log(u) >= 0.5 * x * x + d * (1 - v + log(v)));

    return b * d * v;
}


// Random beta.
double rng_beta(rng_t* rng, double a, double b)
{
    double u = rng_gamma(rng, a, 1.0);
    double v = rng_gamma(rng, b, 1.0);
    return u / (u + v);
}


// Random binomial.
uint32_t rng_binomial(rng_t* rng, double p, uint32_t n)
{
    uint32_t i, a, b, k = 0;
    while (n > 10) {
        a = 1 + (n / 2);
        b = 1 + n - a;

        double x = rng_beta(rng, (double) a, (double) b);
        if (x >= p) {
            n = a - 1;
            p /= x;
        }
        else {
            k += a;
            n = b -1;
            p = (p - x) / (1 - x);
        }
    }

    for (i = 0; i < n; ++i) {
        double u = rng_uniform(rng);
        if (u < p) ++k;
    }

    return k;
}


uint32_t rng_hypergeometric(rng_t* rng, uint32_t n1, uint32_t n2, uint32_t t)
{
    const uint32_t n = n1 + n2;

    uint32_t a = n1;
    uint32_t b = n1 + n2;
    uint32_t k = 0;

    if (t > n) t = n;

    if (t < n/2) {
        for (uint32_t i = 0; i < t; ++i) {
            double u = rng_uniform(rng);
            if (b * u < a) {
                ++k;
                if (k == n1) return k;
                --a;
            }
            --b;
        }
        return k;
    }
    else {
        for (uint32_t i = 0; i < n - t; ++i) {
            double u = rng_uniform(rng);
            if (b * u < a) {
                ++k;
                if (k == n1) return n1 - k;
                --a;
            }
            --b;
        }
        return n1 - k;
    }
}


// Fisher-Yates shuffle
void shuffle(rng_t* rng, uint32_t* ks, uint32_t n)
{
    uint32_t i, j, k;
    for (i = n - 1; i > 0; --i) {
        j = rng_uniform_int(rng, i + 1);
        k = ks[j]; ks[j] = ks[i]; ks[i] = k; // swap
    }
}


// Simple bit vectors
// ------------------

typedef struct bitset_t
{
    // Size in bits.
    size_t n;

    // Size in array elements. (Which is ceiling(n / 64))
    size_t size;

    // Where the bits live.
    uint64_t* xs;
} bitset_t;


void bitset_init(bitset_t* s, size_t n)
{
    s->n = n;
    s->size = (n + 64 - 1) / 64;
    s->xs = malloc_or_die(s->size * sizeof(uint64_t));
    memset(s->xs, 0, s->size * sizeof(uint64_t));
}


void bitset_free(bitset_t* s)
{
    if (s) {
        free(s->xs);
    }
}


void bitset_set(bitset_t* s, size_t i)
{
    s->xs[i / 64] |=  UINT64_C(1) << (i % 64);
}


bool bitset_get(bitset_t* s, size_t i)
{
    return (s->xs[i / 64] >> (i % 64)) & UINT64_C(0x1);
}


// A few notes on sampling without replacement with minimal memory usage,
// assuming both n and m are large. We could use reservoir sampling and
// avoid making two passes over the files, but this would involve allocating
// enough memory to hold m lines. If m is large and the lines are very long,
// this can be prohibitive. Instead, we make one pass to get the total
// number of lines, select samples lines by setting bits it a bitset, then
// make another pass to print the selected lines.
//
// There are multiple strategies we could use to initialize the bitset. The
// easiest is to create an array of indexes 1..n, randomly shuffle them, and
// use the first n indexes to set bits. But this defeats the purpose, since
// the array of indexes would be 8n bytes. We instead adopt a divide and
// conqueror approach. An interval of the bitset is divided in two, we
// allocate n into u + v = n, and recursively set u random bits in the first
// half and v in the second half. If the chunk size is small enough, we use
// the fisher-yates trick. So, this is O(nlog(n)), as opposed to just using
// fisher-yates which is O(n), but if the chunk size is reasonably large,
// performance shouldn't be an issue.

// Resort to the fisher-yates method when the chunk size is this small or less.
const uint64_t brute_chunk_size = 1024;

// Work space for the fisher yates method.
uint32_t ks[brute_chunk_size];

// Set n random bits in the interval [a, b] in the given bitset.
void random_bits(rng_t* rng, bitset_t* s, uint64_t a, uint64_t b, uint64_t n)
{
    assert(a <= b);
    uint64_t m = b - a + 1;
    assert(n <= m);
    if (m <= brute_chunk_size) {
        for (uint32_t i = 0; i < m; ++i) ks[i] = i;
        shuffle(rng, ks, m);
        for (uint32_t i = 0; i < n; ++i) bitset_set(s, a + ks[i]);
    }
    else {
        uint64_t m1 = m / 2;
        uint32_t n1 = rng_hypergeometric(rng, m1, m - m1, n);
        random_bits(rng, s, a, a + m1 - 1, n1);
        random_bits(rng, s, a + m1, b, n - n1);
    }
}


// Count the number of lines in a file.
uint64_t count_lines(FILE* file)
{
    char buf[1024];
    size_t linecount = 0;
    size_t n;
    bool nonempty = false; // is the last line non-empty
    while ((n = fread(buf, 1, sizeof(buf), file))) {
        nonempty = true;
        char* c = memchr(buf, '\n', n);
        while (c != NULL) {
            ++linecount;
            c = memchr(c + 1, '\n', n - (c - buf) - 1);
            nonempty = c && c + 1 < buf + n;
        }
    }
    return nonempty ? linecount + 1 : linecount;
}


void print_help(FILE* out)
{
    fprintf(out, "Usage: subsample [options] in_file [in_file2...] > out_file\n");
}


int main(int argc, char* argv[])
{
    uint32_t seed = UINT32_C(801239084);
    uint32_t n = 1; // return a single random entry by default.
    double p = NAN;

    static struct option long_options[] = {
        {"seed", required_argument, NULL, 's'},
        {"help", no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };

    while (true) {
        int optidx;
        int opt = getopt_long(argc, argv, "n:p:s:h", long_options, &optidx);

        if (opt == -1) break;
        else if (opt == 'n') {
            n = strtoul(optarg, NULL, 10);
        }
        else if (opt == 'p') {
            p = atof(optarg);
            if (p < 0.0) {
                fputs("Sample proportion ('-p') is less than zero.\n", stderr);
                return EXIT_FAILURE;
            }
        }
        else if (opt == 's') {
            if (optarg == NULL) {
                seed = (uint32_t) time(NULL);
            }
            else {
                seed = strtoul(optarg, NULL, 10);
            }
        }
        else if (opt == 'h') {
            print_help(stdout);
            return EXIT_SUCCESS;
        }
        else {
            return EXIT_FAILURE;
        }
    }

    if (optind >= argc) {
        print_help(stderr);
        return EXIT_FAILURE;
    }

    // open all files
    char** filenames = &argv[optind];
    size_t numfiles = argc - optind;
    FILE** files = malloc_or_die(numfiles * sizeof(FILE*));
    for (size_t i = 0; i < numfiles; ++i) {
        files[i] = fopen(filenames[i], "r");
        if (!files[i]) {
            fprintf(stderr, "Error: cannot open %s for reading.\n",
                    filenames[i]);
            return EXIT_FAILURE;
        }
    }

    // count total lines
    uint64_t m = 0;
    for (size_t i = 0; i < numfiles; ++i) {
        m += count_lines(files[i]);
        if (fseek(files[i], 0L, SEEK_SET) != 0) {
            fprintf(stderr, "Error: %s is unseekable.\n", filenames[i]);
            return EXIT_FAILURE;
        }
    }

    if (n > m) {
        fprintf(stderr,
               "Error: cannot sample %lu lines from %lu without replacement\n",
               (unsigned long) n, (unsigned long) m);
        return EXIT_FAILURE;
    }

    rng_t rng;
    rng_init(&rng, seed);

    bitset_t bitset;
    bitset_init(&bitset, m);
    random_bits(&rng, &bitset, 0, m - 1, n);

    char buf[1024];
    size_t line = 0;
    bool sampled = bitset_get(&bitset, 0);
    for (size_t i = 0; i < numfiles; ++i) {
        size_t readlen;
        while ((readlen = fread(buf, 1, sizeof(buf), files[i]))) {
            char* c = buf;
            while (c != NULL && c < buf + readlen) {
                char* cnext = memchr(c, '\n', readlen - (c - buf));
                if (cnext == NULL) {
                    if (sampled) fwrite(c, 1, readlen - (c - buf), stdout);
                    c = NULL;
                }
                else {
                    if (sampled) fwrite(c, 1, cnext - c + 1, stdout);
                    sampled = bitset_get(&bitset, ++line);
                    c = cnext + 1;
                }
            }
        }
    }

    bitset_free(&bitset);
    for (size_t i = 0; i < numfiles; ++i) {
        fclose(files[i]);
    }
}





