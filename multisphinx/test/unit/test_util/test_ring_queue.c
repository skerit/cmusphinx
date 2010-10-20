#include <stdio.h>

#include <sphinxbase/ring_queue.h>

#include "bptbl.h"

int
main(int argc, char *argv[])
{
	ring_queue_t *a, *b;

	a = ring_queue_init(10, sizeof(bp_t), 0);
	b = ring_queue_init(10, sizeof(bp_t), 0);
}
