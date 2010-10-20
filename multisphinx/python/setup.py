try:
    from setuptools import setup, Extension
except:
    from distutils.core import setup, Extension

import distutils.command.install
import os
import commands

class bogus_uninstall(distutils.command.install.install):
    """
    Slightly bogus uninstall, just here to satisfy automake's make
    distcheck.  Do NOT actually use this to uninstall the module!
    """
    def run(self):
        # I believe the word 'bogus' is operative here.  When we run
        # get_outputs() it will create subcommands, which will try to
        # create the original 'install' object, which does not exist
        # at this point.  We need to make sure that the --prefix
        # argument gets propagated to said object.  This is not the
        # right way to do that, but it works, for now.
        install = self.distribution.get_command_obj('install')
        install.prefix = self.prefix
        install.ensure_finalized()
        dirs = {}
        for f in self.get_outputs():
            dirs[os.path.dirname(f)] = 1
            if os.path.isdir(f):
		dirs[f] = 1
		continue
            print "Trying to remove file", f
            try:
                os.unlink(f)
            except:
                pass
        # Gently try to remove any empty directories.
        # This is really not guaranteed to work!!!
        for d in dirs:
            while d != self.prefix:
                print "Trying to remove dir", d
                try:
		    if d.endswith(".egg-info"):
			files=[os.path.join(d,f) for f in os.listdir(d)]
			print "Trying to remove:", " ".join(files)
			for f in files: os.unlink(f)
                    os.rmdir(d)
                except:
                    pass
                d = os.path.dirname(d)

def pkgconfig(*packages, **kw):
    flag_map = {'-I': 'include_dirs', '-L': 'library_dirs', '-l': 'libraries'}
    for token in commands.getoutput("pkg-config --libs --cflags %s" % ' '.join(packages)).split():
        kw.setdefault(flag_map.get(token[:2]), []).append(token[2:])
    return kw

# We actually just need the header files so that we know how to unbox things
try:
    pygtk_includes = pkgconfig('pygtk-2.0')['include_dirs']
except KeyError:
    pygtk_includes = []

setup(name = 'PocketSphinx',
      version = '0.0',
      author = 'David Huggins-Daines',
      author_email = 'dhuggins@cs.cmu.edu',
      description = 'Python interface to CMU PocketSphinx speech recognition',
      license = 'BSD',
      url = 'http://cmusphinx.sourceforge.net',
      ext_modules = [
        Extension('pocketsphinx',
                   sources=['pocketsphinx.c'],
                   libraries=['pocketsphinx', 'sphinxbase'],
                   include_dirs=['../include',
                                 '../include',
                                 '../python',
                                 ] + pygtk_includes,
                   library_dirs=['../src/libpocketsphinx/.libs',
                                 '../src/libsphinxbase/.libs'
                                 ]),
        Extension('sphinxbase',
                   sources=['sphinxbase.c'],
                   libraries=['sphinxbase'],
                   include_dirs=['../include',
                                 '../include'],
                   library_dirs=['../src/libsphinxbase/.libs'])
        ],
      cmdclass = {'bogus_uninstall' : bogus_uninstall}
     ) 