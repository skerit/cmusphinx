SUBDIRS = sphinxbase multisphinx . programs swig test

EXTRA_DIST = autogen.sh \
	m4/pkg.m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = multisphinx.pc
CLEANFILES = multisphinx.pc

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libmultisphinx.la

libmultisphinx_la_SOURCES =
libmultisphinx_la_LDFLAGS = -version-info 0:0:0
libmultisphinx_la_LIBADD = $(LTLIBICONV) \
	sphinxbase/libsphinxbase.la \
	multisphinx/libmultisphinx.la
