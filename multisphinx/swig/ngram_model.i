%feature("docstring","
    N-Gram language model class.

    This class provides access to N-Gram language models stored on
    disk.  These can be in ARPABO text format or Sphinx DMP format.
    Methods are provided for scoring N-Grams based on the model
    and looking up words in the model.

    @param file: Path to an N-Gram model file.
    @type file: string
    @param lw: Language weight to apply to model probabilities.
    @type lw: float
    @param wip: Word insertion penalty to add to model probabilities
    @type wip: float
    @param uw: Weight to give unigrams when interpolating with uniform distribution.
    @type uw: float
") NGramModel;
typedef struct ngram_model_s {
} NGramModel;

%extend NGramModel {
	NGramModel(char const *file, float lw=1.0, float wip=1.0,
		   float uw=1.0) {
		NGramModel *lm;
		logmath_t *lmath = logmath_init(1.0001, 0, 0);

		lm = ngram_model_read(NULL, file, NGRAM_AUTO, lmath);
		logmath_free(lmath);
		if (lw != 1.0 || wip != 1.0 || uw != 1.0)
			ngram_model_apply_weights(lm, lw, wip, uw);
		return lm;
	}
	~NGramModel() {
		ngram_model_free($self);
	}
	%feature("docstring","
        Change the language model weights applied in L{score}.
        
        @param lw: Language weight to apply to model probabilities.
        @type lw: float
        @param wip: Word insertion penalty to add to model probabilities
        @type wip: float
        @param uw: Weight to give unigrams when interpolating with uniform distribution.
        @type uw: float
") apply_weights;
	void apply_weights(float lw=1.0, float wip=1.0, float uw=1.0) {
		ngram_model_apply_weights($self, lw, wip, uw);
	}
	%feature("docstring","
        Get the order of this model (i.e. the 'N' in 'N-gram')

        @return: Order of this model
        @rtype: int
") get_size;
	int get_size() {
		return ngram_model_get_size($self);
	}
};

