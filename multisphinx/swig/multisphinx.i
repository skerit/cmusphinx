%module multisphinx
#ifdef SWIGJAVA
%include "pocketsphinx_java.i"
#endif
#ifdef SWIGPYTHON
%include "argcargv.i"
%include "pocketsphinx_python.i"
#endif
%include "common.i"
%include "carrays.i"

%{
#include <multisphinx/dict.h>
#include <multisphinx/dict2pid.h>
#include <multisphinx/ngram_trie.h>

#include <pocketsphinx.h>

/* Definitions for C's benefit. */
typedef struct logmath_s LogMath;
typedef struct cmd_ln_s Config;
typedef struct bin_mdef_s Mdef;
typedef struct ngram_model_s NGramModel;
typedef struct ngram_iter_s NGramIterator;
typedef dict_t Dict;
typedef dict2pid_t DictToPid;
typedef struct ngram_trie_s NGramTrie;
typedef xwdssid_t xwdssid;

%}


/* Now define them again for SWIG. */
%feature("docstring", "
    Log-space math class.
    
    This class provides fast logarithmic math functions in base
    1.000+epsilon, useful for fixed point speech recognition.

    @param base: The base B in which computation is to be done.
    @type base: float
    @param shift: Log values are shifted right by this many bits.
    @type shift: int
    @param use_table Whether to use an add table or not
    @type use_table: bool") LogMath;
typedef struct logmath_s {
} LogMath;

%extend LogMath {
	LogMath(double base, int shift=0, bool use_table=false) {
		return logmath_init(base, shift, use_table);
	}
	~LogMath() {
		logmath_free($self);
	}
	%feature("docstring","
        Get the log base.

        @return: Logarithmic base used in this object.
        @rtype: float
") get_base;
	double get_base() {
		return logmath_get_base($self);
	}
	%feature("docstring","
        Get the log-zero value.

        @return: Smallest number representable by this object.
        @rtype: int
") get_zero;
	int get_zero() {
		return logmath_get_zero($self);
	}
	%feature("docstring","
        Add two numbers in log-space without using an add-table.

        @param a: Logarithm A.
        @type a: int
        @param b: Logarithm B.
        @type b: int
        @return: log(exp(a)+exp(b))
        @rtype: int
") add;
	int add_exact(int p, int q) {
		return logmath_add_exact($self, p, q);
	}
	%feature("docstring","
        Add two numbers in log-space.

        @param a: Logarithm A.
        @type a: int
        @param b: Logarithm B.
        @type b: int
        @return: log(exp(a)+exp(b))
        @rtype: int
") add;
	int add(int p, int q) {
		return logmath_add($self, p, q);
	}
	%feature("docstring","
        Return log-value of a number.

        @param x: Number (in linear space)
        @type x: float
        @return: Log-value of x.
        @rtype: int
") log;
	int log(double p) {
		return logmath_log($self, p);
	}
	%feature("docstring", "
        Return linear of a log-value

        @param x: Logarithm X (in this object's base)
        @type x: int
        @return: Exponent (linear value) of X.
        @rtype: float
") exp;
	double exp(int logb_p) {
		return logmath_exp($self, logb_p);
	}
	%feature("docstring", "
        Return log-value of a natural logarithm.

        @param x: Logarithm X (in base e)
        @type x: float
        @return: Log-value equivalent of x.
        @rtype: int
") ln_to_log;
	int ln_to_log(double log_p) {
		return logmath_ln_to_log($self, log_p);
	}
	%feature("docstring", "
        Return natural logarithm of a log-value.

        @param x: Logarithm X (in this object's base)
        @type x: int
        @return: Natural log equivalent of x.
        @rtype: float
") log_to_ln;
	double log_to_ln(int logb_p) {
		return logmath_log_to_ln($self, logb_p);
	}
	%feature("docstring","
        Return log-value of a base 10 logarithm.

        @param x: Logarithm X (in base 10)
        @type x: float
        @return: Log-value equivalent of x.
        @rtype: int
") log10_to_log;
	int log10_to_log(double log_p) {
		return logmath_log10_to_log($self, log_p);
	}
	%feature("docstring","
        Return logarithm in base 10 of a log-value.

        @param x: Logarithm X (in this object's base)
        @type x: int
        @return: log10 equivalent of x.
        @rtype: float
") log_to_log10;
	double log_to_log10(int logb_p) {
		return logmath_log_to_log10($self, logb_p);
	}
}

typedef struct cmd_ln_s {
} Config;

%extend Config {
	Config() {
		Config *c = cmd_ln_init(NULL, ps_args(), FALSE, NULL);
		return c;
	}
	Config(char const *file) {
		Config *c = cmd_ln_parse_file_r(NULL, ps_args(), file, FALSE);
		return c;
	}
	Config(int ARGC, char **ARGV) {
		Config *c = cmd_ln_parse_r(NULL, ps_args(), ARGC, ARGV, FALSE);
		return c;
	}
	~Config() {
		cmd_ln_free_r($self);
	}
	void setBoolean(char const *key, bool val) {
		cmd_ln_set_boolean_r($self, key, val);
	}
	void setInt(char const *key, int val) {
		cmd_ln_set_int_r($self, key, val);
	}
	void setFloat(char const *key, double val) {
		cmd_ln_set_float_r($self, key, val);
	}
	void setString(char const *key, char const *val) {
		cmd_ln_set_str_r($self, key, val);
	}
	bool exists(char const *key) {
		return cmd_ln_exists_r($self, key);
	}
	bool getBoolean(char const *key) {
		return cmd_ln_boolean_r($self, key);
	}
	int getInt(char const *key) {
		return cmd_ln_int_r($self, key);
	}
	double getFloat(char const *key) {
		return cmd_ln_float_r($self, key);
	}
	char const *getString(char const *key) {
		return cmd_ln_str_r($self, key);
	}
};

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

typedef struct bin_mdef_s {
} Mdef;

%extend Mdef {
	Mdef(Config *config, char const *filename) {
		return bin_mdef_read(config, filename);
	}
	~Mdef() {
		bin_mdef_free($self);
	}
	int ciphone_id(char const *ciphone) {
		return bin_mdef_ciphone_id($self, ciphone);
	}
	char const *ciphone_str(int ci) {
		return bin_mdef_ciphone_str($self, ci);
	}
	int phone_id(int b, int l, int r, char pos) {
		return bin_mdef_phone_id($self, b, l, r, pos);
	}
	int phone_id_nearest(int b, int l, int r, char pos) {
		return bin_mdef_phone_id_nearest($self, b, l, r, pos);
	}
	bool is_fillerphone(int pid) {
		return bin_mdef_is_fillerphone($self, pid);
	}
	bool is_ciphone(int pid) {
		return bin_mdef_is_ciphone($self, pid);
	}
	int n_ciphone() {
		return bin_mdef_n_ciphone($self);
	}
	int n_phone() {
		return bin_mdef_n_phone($self);
	}
	int n_sseq() {
		return bin_mdef_n_sseq($self);
	}
	int silphone() {
		return bin_mdef_silphone($self);
	}
}

typedef struct {
} Dict;

%extend Dict {
	Dict(Mdef *mdef, char const *dictfile=NULL, char const *fdictfile=NULL) {
		Config *config = cmd_ln_init(NULL, ps_args(), FALSE, NULL);
		Dict *d;
		if (dictfile) cmd_ln_set_str_r(config, "-dict", dictfile);
		if (fdictfile) cmd_ln_set_str_r(config, "-fdict", fdictfile);
		d = dict_init(config, mdef);
		cmd_ln_free_r(config);
		return d;
	}
	~Dict() {
		dict_free($self);
	}
	int wordid(char const *word) {
		return dict_wordid($self, word);
	}
	bool filler_word(int wid) {
		return dict_filler_word($self, wid);
	}
	bool real_word(int wid) {
		return dict_real_word($self, wid);
	}
	char const *wordstr(int wid) {
		return dict_wordstr($self, wid);
	}
	char const *basestr(int wid) {
		return dict_basestr($self, wid);
	}
	int basewid(int wid) {
		return dict_basewid($self, wid);
	}
	int pronlen(int wid) {
		return dict_pronlen($self, wid);
	}
	int first_phone(int wid) {
		return dict_first_phone($self, wid);
	}
	int second_phone(int wid) {
		return dict_second_phone($self, wid);
	}
	int second_last_phone(int wid) {
		return dict_second_last_phone($self, wid);
	}
	int last_phone(int wid) {
		return dict_last_phone($self, wid);
	}
	int size() {
		return dict_size($self);
	}
	int num_fillers() {
		return dict_num_fillers($self);
	}
	int startwid() {
		return dict_startwid($self);
	}
	int finishwid() {
		return dict_finishwid($self);
	}
	int silwid() {
		return dict_silwid($self);
	}
	/* FIXME: Pythonify this, eventually. */
	int pron(int wid, int pos) {
		return dict_pron($self, wid, pos);
	}
}

typedef int s3ssid_t;
typedef int s3cipid_t;

%array_class(s3ssid_t, ssidArray)
%array_class(s3cipid_t, cipidArray)

typedef struct {
	ssidArray *ssid;
	cipidArray *cimap;
	int n_ssid;
} xwdssid;

typedef struct {
} DictToPid;

%extend DictToPid {
	DictToPid(Mdef *mdef, Dict *dict) {
		return dict2pid_build(mdef, dict);
	}
	~DictToPid() {
		dict2pid_free($self);
	}
	Mdef* get_mdef() {
		return bin_mdef_retain($self->mdef);
	}
	Dict* get_dict() {
		return dict_retain($self->dict);
	}
	xwdssid *rssid(int ci, int lc) {
		return dict2pid_rssid($self, ci, lc);
	}
}

typedef struct ngram_trie_s {
} NGramTrie;

%extend NGramTrie {
	NGramTrie(Dict *d, LogMath *lmath) {
		return ngram_trie_init(d, lmath);
	}
	~NGramTrie() {
		ngram_trie_free($self);
	}
}
