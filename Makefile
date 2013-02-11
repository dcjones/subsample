
CFLAGS=-Wall -Wextra -g -O3 -std=c99

subsample: subsample.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f subsample


