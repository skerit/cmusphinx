
typedef struct ps_search_s {
} Search;

%extend Search {
	Search() {
		return NULL;
	}
	~Search() {
		ps_search_free($self);
	}
	void run() {
		ps_search_run($self);
	}
	void wait() {
		ps_search_wait($self);
	}
	char const *hyp(int *OUTPUT) {
		return ps_search_hyp($self, OUTPUT);
	}
};

Search *fwdtree_search_init(Config *c, Acmod *am, Dict *d, DictToPid *d2p);
Search *fwdflat_search_init(Config *c, Acmod *am, Dict *d, DictToPid *d2p,
			    NGramModel *lmset=NULL);
/* FIXME (if possible??!): There is no type checking here so these
 * functions will segfault if you look at them funny. */
VocabMap *fwdflat_search_set_vocab_map(Search *s, VocabMap *vm);

