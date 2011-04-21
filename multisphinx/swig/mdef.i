
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
};

