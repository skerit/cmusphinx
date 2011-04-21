
typedef struct bptbl_s {
} Bptbl;

%extend Bptbl {
	Bptbl(char const *name, DictToPid *d2p,
	      int n_alloc=512, int n_frame_alloc=128) {
		return bptbl_init(name, d2p, n_alloc, n_frame_alloc);
	}
	~Bptbl() {
		bptbl_free($self);
	}
};

