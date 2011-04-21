typedef struct vocab_map_s {
} VocabMap;

%extend VocabMap {
	VocabMap(Dict *d) {
		return vocab_map_init(d);
	}
	~VocabMap() {
		vocab_map_free($self);
	}
};
