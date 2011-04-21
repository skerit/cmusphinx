SUBDIRS = \
	libsphinxbase \
	libpocketsphinx

lib_LTLIBRARIES = libmultisphinx.la
libmultisphinx_la_LDFLAGS = -version-info 0:0:0
libmultisphinx_la_SOURCES =
libmultisphinx_la_LIBADD = $(LTLIBICONV) \
	libsphinxbase/util/libsphinxutil.la \
	libsphinxbase/fe/libsphinxfe.la \
	libsphinxbase/feat/libsphinxfeat.la \
	libsphinxbase/lm/libsphinxlm.la \
	libpocketsphinx/libpocketsphinx.la
