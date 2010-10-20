SUBDIRS = util fe feat lm

lib_LTLIBRARIES = libsphinxbase.la

libsphinxbase_la_SOURCES =
# Confusing libtool version info means:
#
# This is API version 3, software release 0, which only supports API 3
libsphinxbase_la_LDFLAGS = -version-info 3;0:0
libsphinxbase_la_LIBADD = $(LTLIBICONV) \
	util/libsphinxutil.la \
	fe/libsphinxfe.la \
	feat/libsphinxfeat.la \
	lm/libsphinxlm.la
