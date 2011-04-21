
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
};
