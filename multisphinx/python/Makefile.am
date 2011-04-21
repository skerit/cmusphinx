EXTRA_DIST = pocketsphinx.c \
	pocketsphinx.pyx \
	pocketsphinx.pxd \
	sphinxbase.c \
	sphinxbase.pyx \
	sphinxbase.pxd \
	sb_test.py \
	hufftest2.py \
	ps_test.py \
	setup_win32.py

pkginclude_HEADERS = pocketsphinx.pxd sphinxbase.pxd
noinst_HEADERS = bogus_pygobject.h

if BUILD_CYTHON
$(srcdir)/pocketsphinx.c: $(srcdir)/pocketsphinx.pyx $(srcdir)/pocketsphinx.pxd
	cython -o $@ $< $(CPPFLAGS) -I$(sphinxbase)/python
$(srcdir)/sphinxbase.c: $(srcdir)/sphinxbase.pyx $(srcdir)/sphinxbase.pxd
	cython -o $@ $< $(CPPFLAGS) -I$(sphinxbase)/python
endif
