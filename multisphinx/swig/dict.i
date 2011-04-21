%{
#include <multisphinx/cmdln_macro.h>
static arg_t const dict_args[] = {
	DICT_OPTIONS,
	CMDLN_EMPTY_OPTION
};
%}

typedef struct {
} Dict;

%extend Dict {
	Dict(Mdef *mdef, char const *dictfile=NULL, char const *fdictfile=NULL) {
		cmd_ln_t *config = cmd_ln_init(NULL, dict_args, FALSE, NULL);
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
};
