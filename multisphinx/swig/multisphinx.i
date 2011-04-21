%module multisphinx
#ifdef SWIGJAVA
%include "includes_java.i"
#endif
#ifdef SWIGPYTHON
%include "argcargv.i"
%include "includes_python.i"
#endif
%include "common.i"
%include "carrays.i"
%include "cdata.i"

%{
#include <multisphinx/dict.h>
#include <multisphinx/dict2pid.h>
#include <multisphinx/ngram_trie.h>
#include <multisphinx/featbuf.h>
#include <multisphinx/bptbl.h>
#include <multisphinx/arc_buffer.h>
#include <multisphinx/fwdtree_search.h>
#include <multisphinx/fwdflat_search.h>
#include <multisphinx/pocketsphinx.h>

/* Definitions for C's benefit. */
typedef struct logmath_s LogMath;
typedef struct cmd_ln_s Config;
typedef struct bin_mdef_s Mdef;
typedef struct ngram_model_s NGramModel;
typedef struct ngram_iter_s NGramIterator;
typedef struct featbuf_s FeatBuf;
typedef struct ps_search_s Search;
typedef struct acmod_s Acmod;
typedef struct bptbl_s Bptbl;
typedef struct arc_buffer_s ArcBuffer;
typedef struct vocab_map_s VocabMap;
typedef dict_t Dict;
typedef dict2pid_t DictToPid;
typedef struct ngram_trie_s NGramTrie;
typedef xwdssid_t xwdssid;

%}

/* Now define them again for SWIG. */
%include "logmath.i"
%include "cmd_ln.i"
%include "ngram_model.i"
%include "mdef.i"
%include "dict.i"
%include "dict2pid.i"
%include "ngram_trie.i"
%include "featbuf.i"
%include "acmod.i"
%include "ps_search.i"
%include "bptbl.i"
%include "arc_buffer.i"
%include "vocab_map.i"
