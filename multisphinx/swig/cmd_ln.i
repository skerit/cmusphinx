%{
#include <multisphinx/cmdln_macro.h>
static arg_t const ps_args[] = {
	MULTISPHINX_OPTIONS,
	CMDLN_EMPTY_OPTION
};
%}

typedef struct cmd_ln_s {
} Config;

%extend Config {
	Config() {
		Config *c = cmd_ln_init(NULL, ps_args, FALSE, NULL);
		return c;
	}
	Config(char const *file) {
		Config *c = cmd_ln_parse_file_r(NULL, ps_args, file, FALSE);
		return c;
	}
	Config(int ARGC, char **ARGV) {
		Config *c = cmd_ln_parse_r(NULL, ps_args, ARGC, ARGV, FALSE);
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
