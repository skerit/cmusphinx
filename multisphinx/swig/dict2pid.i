
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
};

