
typedef struct arc_buffer_s {
} ArcBuffer;

%extend ArcBuffer {
	ArcBuffer(char const *name, Bptbl *input_bptbl,
		  NGramModel *lm, bool keep_scores) {
		return arc_buffer_init(name, input_bptbl, lm, keep_scores);
	}
	~ArcBuffer() {
		arc_buffer_free($self);
	}
};

