#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <multisphinx/ms_lattice.h>

#include "test_macros.h"

int
main(int argc, char *argv[])
{
	ms_lattice_t *l;
	ms_latnode_t *s, *a, *b, *c, *d, *e;
	ms_latnode_iter_t *itor;
	logmath_t *lmath;
	dict_t *dict;
	FILE *fh;
	int32 idx, idx2;

	lmath = logmath_init(1.0001, 0, FALSE);
	TEST_ASSERT(l = ms_lattice_init(lmath, NULL));
	dict = ms_lattice_dict(l);

	/* Add some nodes and arcs and stuff. */
	idx = ms_lattice_lmstate_init(l, dict_startwid(dict), NULL, 0);
	printf("idx(<s>)=%d\n", idx);
	s = ms_lattice_node_init(l, 0, idx);
	ms_lattice_set_start(l, s);
	TEST_ASSERT(s == ms_lattice_get_node_id
		    (l, 0,
		     ms_lattice_get_lmstate_idx(l, dict_startwid(dict), NULL, 0)));

	/* <s> -> A : <s>/-42 */
	idx2 = ms_lattice_lmstate_init(l, dict_add_word(dict, "A", NULL, 0),
				       NULL, 0);
	printf("idx(A)=%d\n", idx2);
	a = ms_lattice_node_init(l, 5, idx2);
	ms_lattice_link(l, s, a, idx, -42);
	TEST_ASSERT(a == ms_lattice_get_node_id
		    (l, 5,
		     ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "A"), NULL, 0)));

	/* <s> -> B : <s>/-69 */
	idx2 = ms_lattice_lmstate_init(l, dict_add_word(dict, "B", NULL, 0),
				       NULL, 0);
	printf("idx(B)=%d\n", idx2);
	b = ms_lattice_node_init(l, 7, idx2);
	ms_lattice_link(l, s, b, idx, -69);
	TEST_ASSERT(b == ms_lattice_get_node_id
		    (l, 7,
		     ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "B"), NULL, 0)));

	/* B -> C : B/-420 */
	idx = ms_lattice_lmstate_init(l, dict_add_word(dict, "C", NULL, 0), NULL, 0);
	printf("idx(C)=%d\n", idx);
	c = ms_lattice_node_init(l, 8, idx);
	ms_lattice_link(l, b, c, idx2, -420);
	TEST_ASSERT(c == ms_lattice_get_node_id
		    (l, 8,
		     ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "C"), NULL, 0)));

	/* B -> D : B/-420 */
	idx = ms_lattice_lmstate_init(l, dict_add_word(dict, "D", NULL, 0), NULL, 0);
	printf("idx(D)=%d\n", idx);
	d = ms_lattice_node_init(l, 9, idx);
	ms_lattice_link(l, b, d, idx2, -666);
	TEST_ASSERT(d == ms_lattice_get_node_id
		    (l, 9,
		     ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "D"), NULL, 0)));

	/* A -> D : A/-420 */
	idx2 = ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "A"), NULL, 0);
	printf("idx()=%d\n", idx2);
	ms_lattice_link(l, a, d, idx2, -999);

	/* C -> </s> : C/-99 */
	idx = ms_lattice_lmstate_init(l, dict_finishwid(dict), NULL, 0);
	printf("idx(</s>)=%d\n", idx);
	e = ms_lattice_node_init(l, 13, idx);
	TEST_ASSERT(e == ms_lattice_get_node_id
		    (l, 13,
		     ms_lattice_get_lmstate_idx(l, dict_finishwid(dict), NULL, 0)));
	idx2 = ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "C"), NULL, 0);
	printf("idx(C)=%d\n", idx2);
	ms_lattice_link(l, c, e, idx2, -99);

	/* D -> </s> : D/-69 */
	idx2 = ms_lattice_get_lmstate_idx(l, dict_wordid(dict, "D"), NULL, 0);
	printf("idx(D)=%d\n", idx2);
	ms_lattice_link(l, d, e, idx2, -69);
	ms_lattice_set_end(l, e);

	/* Write it out to a dot file so we can verify it. */
	TEST_ASSERT(fh = fopen("test_build_lattice.dot", "w"));
	TEST_ASSERT(0 == ms_lattice_write_dot(l, fh));
	TEST_ASSERT(0 == fclose(fh))

	/* Test traversal functions. */
	itor = ms_lattice_traverse_topo(l, NULL);
	TEST_ASSERT(ms_latnode_iter_get(itor) == s);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == a);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == b);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == c);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == d);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == e);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_traverse_topo(l, d);
	TEST_ASSERT(ms_latnode_iter_get(itor) == s);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == a);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == b);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == c);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == d);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_reverse_topo(l, NULL);
	TEST_ASSERT(ms_latnode_iter_get(itor) == e);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == c);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == d);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == b);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == a);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == s);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_reverse_topo(l, b);
	TEST_ASSERT(ms_latnode_iter_get(itor) == e);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == c);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == d);
	TEST_ASSERT(ms_latnode_iter_next(itor));
	TEST_ASSERT(ms_latnode_iter_get(itor) == b);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_traverse_frame(l, 0);
	TEST_ASSERT(ms_latnode_iter_get(itor) == s);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_traverse_frame(l, 8);
	TEST_ASSERT(ms_latnode_iter_get(itor) == c);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_traverse_frame(l, 13);
	TEST_ASSERT(ms_latnode_iter_get(itor) == e);
	TEST_ASSERT(ms_latnode_iter_next(itor) == NULL);

	itor = ms_lattice_traverse_frame(l, 2);
	TEST_ASSERT(itor == NULL);

	itor = ms_lattice_traverse_frame(l, 999);
	TEST_ASSERT(itor == NULL);

	ms_lattice_free(l);
	return 0;
}
