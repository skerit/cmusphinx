headers = $(top_srcdir)/include/pocketsphinx.h \
	$(top_srcdir)/include/fsg_set.h \
	$(top_srcdir)/include/ps_lattice.h

latex/refman.pdf: doxyfile $(headers)
	doxygen
	$(MAKE) -C latex refman.pdf

html/index.html: doxyfile $(headers)
	doxygen

clean-local:
	-rm -rf html latex doxytags
