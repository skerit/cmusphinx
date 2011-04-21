bin_PROGRAMS = \
	multisphinx_batch \
	multisphinx_map_vocab \
	multisphinx_merge_lmclasses

LDADD = \
	$(top_builddir)/src/libsphinxbase/libsphinxbase.la \
	$(top_builddir)/src/libpocketsphinx/libpocketsphinx.la
INCLUDES = -I$(top_srcdir)/include \
	-I$(top_srcdir)/include/multisphinx \
	-I$(top_srcdir)/src/libpocketsphinx \
	-I$(top_builddir)/include
