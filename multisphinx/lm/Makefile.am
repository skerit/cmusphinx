noinst_LTLIBRARIES = libsphinxlm.la

AM_YFLAGS = -d
BUILT_SOURCES = jsgf_parser.h

libsphinxlm_la_SOURCES =			\
	lm3g_model.c				\
	ngram_model.c				\
	ngram_model_arpa.c			\
	ngram_model_dmp.c			\
	ngram_model_dmp32.c			\
	ngram_model_set.c			\
	fsg_model.c				\
	jsgf.c					\
	jsgf_scanner.l				\
	jsgf_parser.y

libsphinxlm_la_LIBADD = $(LTLIBICONV)

noinst_HEADERS = ngram_model_internal.h		\
	ngram_model_dmp.h			\
	ngram_model_set.h			\
	ngram_model_arpa.h			\
	lm3g_model.h				\
	jsgf_internal.h				\
	jsgf_scanner.h				\
	jsgf_parser.h

INCLUDES = -I$(top_srcdir)/include \
	   -I$(top_srcdir)/include/sphinxbase \
           -I$(top_builddir)/include 

EXTRA_DIST = lm3g_templates.c
