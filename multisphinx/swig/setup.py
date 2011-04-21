from setuptools import setup, Extension

multisphinx = Extension("_multisphinx",
                        sources = [ "multisphinx.i" ],
                        libraries = [ "multisphinx" ],
                        include_dirs = [".."],
                        library_dirs = ["../.libs"],
                        depends = [ "acmod.i", "arc_buffer.i", "bptbl.i", "cmd_ln.i",
                                   "common.i", "dict.i", "dict2pid.i", "featbuf.i",
                                   "includes_python.i", "logmath.i", "mdef.i",
                                   "ngram_model.i", "ngram_trie.i", "search_factory.i",
                                   "search.i", "vocab_map.i", "state_align.i"])

setup(name = "multisphinx",
      version = "0.2",
      ext_modules = [multisphinx],
      py_modules = ["multisphinx"])
