typedef struct {
	union {
        int wid;
        struct {
            int ssid;
            int cipid;
            int tmatid;
        } pid;
        int senid;
	} id;
	int start;
	int duration;
} AlignmentEntry;

typedef struct {
} AlignmentIter;

%nodefaultctor AlignmentIter;
%extend AlignmentIter {
	~AlignmentIter() {
		alignment_iter_free($self);
	}
	/* next() is special for Python, Java, etc... */
	AlignmentIter *_next() {
		return alignment_iter_next($self);
	}
	AlignmentIter *prev() {
		return alignment_iter_prev($self);
	}
	AlignmentIter *up() {
		return alignment_iter_up($self);
	}
	AlignmentIter *down() {
		return alignment_iter_down($self);
	}
	AlignmentIter *goto(int idx) {
		return alignment_iter_goto($self, idx);
	}
	AlignmentEntry *get() {
		return alignment_iter_get($self);
	}
%pythoncode %{
def __iter__(self):
	self.starting = True
	self.done = False
	return self

def next(self):
	if self.starting == True:
		self.starting = False
		return self.get()
	if self.done == True:
		spam = None
	else:
		spam = self._next()
	if spam == None:
	    self.done = True
        raise StopIteration
	return self.get() 
%}
}


typedef struct {
} Alignment;

%extend Alignment {
	Alignment(DictToPid *d2p) {
		return alignment_init(d2p);
	}
	~Alignment() {
		alignment_free($self);
	}
	int add_word(char *word, int duration=0) {
		return alignment_add_word($self, dict_wordid($self->d2p->dict, word), duration);
	}
	int add_words(int ARGC, char **ARGV) {
		int i;
		for (i = 0; i < ARGC; ++i) {
			alignment_add_word($self, dict_wordid($self->d2p->dict, ARGV[i]), 0);
		}
		return i;
	}
	int populate() {
		return alignment_populate($self);
	}
	int populate_ci() {
		return alignment_populate_ci($self);
	}
	int propagate() {
		return alignment_propagate($self);
	}
	int n_words() {
		return alignment_n_words($self);
	}
	int n_phones() {
		return alignment_n_phones($self);
	}
	int n_states() {
		return alignment_n_states($self);
	}
	AlignmentIter *words() {
		return alignment_words($self);
	}
	AlignmentIter *phones() {
		return alignment_phones($self);
	}
	AlignmentIter *states() {
		return alignment_states($self);
	}
}

%extend Search {
	int set_alignment(Alignment *al) {
		return state_align_search_set_alignment($self, al);
	}
}