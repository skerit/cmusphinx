
typedef struct acmod_s {
} Acmod;

%extend Acmod {
	Acmod(Config *c, LogMath *lmath, FeatBuf *fb) {
		return acmod_init(c, lmath, fb);
	}
	~Acmod() {
		acmod_free($self);
	}
};
