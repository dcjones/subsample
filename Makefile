
CFLAGS=-Wall -Wextra -g -O0 -std=c99

subsample: subsample.c
	$(CC) $(CFLAGS) -o $@ $^ -lm

clean:
	rm -f subsample


