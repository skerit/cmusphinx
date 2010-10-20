bin_PROGRAMS = pocketsphinx_batch

pocketsphinx_batch_SOURCES = batch.c
pocketsphinx_batch_LDADD = \
	$(top_builddir)/src/libsphinxbase/libsphinxbase.la \
	$(top_builddir)/src/libpocketsphinx/libpocketsphinx.la

INCLUDES = -I$(top_srcdir)/include \
	-I$(top_srcdir)/src/libpocketsphinx \
	-I$(top_builddir)/include
