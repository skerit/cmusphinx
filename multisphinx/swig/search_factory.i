typedef struct {
} SearchFactory;

%extend SearchFactory {
	SearchFactory(int ARGC, char **ARGV) {
		return search_factory_init_argv(ARGC, ARGV);
	}
	~SearchFactory() {
		search_factory_free($self);
	}

	Search *create(char *name) {
		return search_factory_create_argv($self, NULL, name, 0, NULL);
	}

	Search *create(Search *other, char *name) {
		return search_factory_create_argv($self, other, name, 0, NULL);
	}

	Search *create(Search *other, char *name, int ARGC, char **ARGV) {
		return search_factory_create_argv($self, other, name, ARGC, ARGV);
	}

	FeatBuf *featbuf() {
		return search_factory_featbuf($self);
	}

	Acmod *acmod() {
		return search_factory_acmod($self);
	}

	NGramModel *lm() {
		return search_factory_lm($self);
	}

	DictToPid *d2p() {
		return search_factory_d2p($self);
	}

	Dict *dict() {
		return search_factory_d2p($self)->dict;
	}

	Mdef *mdef() {
		return search_factory_d2p($self)->mdef;
	}
}
