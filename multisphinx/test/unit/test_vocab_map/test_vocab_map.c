/**
 * @file test_vocab_map.c
 * @brief Test the vocabulary mapping object.
 * @author David Huggins-Daines <dhuggins@cs.cmu.edu>
 */

#include <pocketsphinx.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "pocketsphinx_internal.h"
#include "vocab_map.h"
#include "test_macros.h"

int
main(int argc, char *argv[])
{
	vocab_map_t *vm;

	vm = vocab_map_init();

	vocab_map_expand_active();
}
