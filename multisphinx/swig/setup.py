from setuptools import setup, Extension

multisphinx = Extension("_multisphinx",
                        sources = [ "multisphinx.i" ],
                        libraries = [ "multisphinx" ],
                        include_dirs = [".."],
                        library_dirs = ["../.libs"])

setup(name = "multisphinx",
      version = "0.1",
      ext_modules = [multisphinx],
      py_modules = ["multisphinx"])
