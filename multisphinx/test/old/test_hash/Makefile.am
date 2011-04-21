check_PROGRAMS = displayhash deletehash test_hash_iter

noinst_HEADERS = test_macros.h

INCLUDES = -I$(top_srcdir)

LDADD = ${top_builddir}/libmultisphinx.la

TESTS = test_hash_iter				\
	_hash_delete1.test			\
	_hash_delete2.test			\
	_hash_delete3.test			\
	_hash_delete4.test			\
	_hash_delete5.test

EXTRA_DIST = display.res			\
	_hash_delete1.res			\
	_hash_delete1.test			\
	_hash_delete2.res			\
	_hash_delete2.test			\
	_hash_delete3.res			\
	_hash_delete3.test			\
	_hash_delete4.res			\
	_hash_delete4.test			\
	_hash_delete5.res			\
	_hash_delete5.test

CLEANFILES = *.out