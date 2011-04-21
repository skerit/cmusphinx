%module pocketsphinx
#ifdef SWIGJAVA
%include "pocketsphinx_java.i"
#endif
#ifdef SWIGPYTHON
%include "pocketsphinx_python.i"
#endif
%include "common.i"

%import "multisphinx.i"

%{
#include <pocketsphinx.h>
#include <sphinxbase/err.h>
#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/cmd_ln.h>

typedef struct cmd_ln_s Config;
typedef struct ps_decoder_s Decoder;

/* Auxiliary objects used to return multiple values. */
typedef struct hyp_s {
	char *hypstr;
	char *uttid;
	int best_score;
} Hypothesis;
%}

/* Auxiliary types for returning multiple values. */
typedef struct hyp_s {
	char *hypstr;
	char *uttid;
	int best_score;
} Hypothesis;
/* These are opaque types but we have to "define" them for SWIG. */
typedef struct ps_decoder_s {
} Decoder;


%extend Hypothesis {
	Hypothesis(char const *hypstr, char const *uttid, int best_score) {
		Hypothesis *h = ckd_calloc(1, sizeof(*h));
		if (hypstr)
			h->hypstr = ckd_salloc(hypstr);
		if (uttid)
			h->uttid = ckd_salloc(uttid);
		h->best_score = best_score;
		return h;
		
	}
	~Hypothesis() {
		ckd_free($self->hypstr);
		ckd_free($self->uttid);
		ckd_free($self);
	}
}

%extend Decoder {
	Decoder() {
		Decoder *d = ps_init(cmd_ln_init(NULL, ps_args(), FALSE, NULL));
		return d;
	}
	Decoder(Config *c) {
		Decoder *d = ps_init(c);
		return d;
	}
	Config *getConfig() {
		return cmd_ln_retain(ps_get_config($self));
	}
	int startUtt() {
		return ps_start_utt($self, NULL);
	}
	int startUtt(char const *uttid) {
		return ps_start_utt($self, uttid);
	}
	char const *getUttid() {
		return ps_get_uttid($self);
	}
	int endUtt() {
		return ps_end_utt($self);
	}
	int processRaw(short const *SDATA, size_t NSAMP, bool no_search, bool full_utt) {
		return ps_process_raw($self, SDATA, NSAMP, no_search, full_utt);
	}
	Hypothesis *getHyp() {
		char const *hyp, *uttid;
		int32 best_score;
		hyp = ps_get_hyp($self, &best_score, &uttid);
		return new_Hypothesis(hyp, uttid, best_score);
	}
	~Decoder() {
		ps_free($self);
	}
};

%inline {
	/* Static method to set the logging file. */
	void setLogfile(char const *path) {
		err_set_logfile(path);
	}
};

