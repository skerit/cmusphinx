bin_PROGRAMS = pocketsphinx_batch \
	multisphinx_batch \
	pocketsphinx_map_vocab \
	ngram_homos

LDADD = \
	$(top_builddir)/src/libsphinxbase/libsphinxbase.la \
	$(top_builddir)/src/libpocketsphinx/libpocketsphinx.la
INCLUDES = -I$(top_srcdir)/include \
	-I$(top_srcdir)/src/libpocketsphinx \
	-I$(top_builddir)/include
