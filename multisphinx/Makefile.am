SUBDIRS = src \
	  doc \
	  include \
	  python

EXTRA_DIST = autogen.sh							\
	swig/Makefile							\
	swig/Android.mk							\
	swig/pocketsphinx.i						\
	m4/pkg.m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = multisphinx.pc
CLEANFILES = multisphinx.pc

ACLOCAL_AMFLAGS = -I m4