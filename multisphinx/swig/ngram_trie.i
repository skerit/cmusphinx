typedef struct ngram_trie_s {
} NGramTrie;

%extend NGramTrie {
	NGramTrie(Dict *d, LogMath *lmath) {
		return ngram_trie_init(d, lmath);
	}
	~NGramTrie() {
		ngram_trie_free($self);
	}
};
