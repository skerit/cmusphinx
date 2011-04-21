from setuptools import setup, Extension

multisphinx = Extension("_multisphinx",
                        sources = [ "multisphinx.i" ],
                        libraries = [ "sphinxbase", "pocketsphinx" ],
                        include_dirs = ["../include"],
                        library_dirs = ["../src/libpocketsphinx/.libs",
                                        "../src/libsphinxbase/.libs"])

setup(name = "multisphinx",
      version = "0.1",
      ext_modules = [multisphinx],
      py_modules = ["multisphinx"])

