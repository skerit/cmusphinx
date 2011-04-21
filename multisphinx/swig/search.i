typedef struct seg_iter_s {
} SegIter;

%extend SegIter {
	SegIter() {
		return NULL;
	}
	~SegIter() {
		seg_iter_free($self);
	}
	SegIter *next() {
		return seg_iter_next($self);
	}
	char *word() {
		return seg_iter_word($self);
	}
	int times(int *OUTPUT) {
		int sf;
		seg_iter_times($self, &sf, OUTPUT);
		return sf;
	}
}

typedef struct search_s {
} Search;

%extend Search {
	Search() {
		return NULL;
	}
	~Search() {
		search_free($self);
	}
	void run() {
		search_run($self);
	}
	void wait() {
		search_wait($self);
	}
	char const *hyp(int *OUTPUT) {
		return search_hyp($self, OUTPUT);
	}
	SegIter *seg_iter() {
		int32 foo;
		return search_seg_iter($self, &foo);
	}
	ArcBuffer *link(Search *other, char *name, bool keep_scores=true) {
		return search_link($self, other, name, keep_scores);
	}
};
