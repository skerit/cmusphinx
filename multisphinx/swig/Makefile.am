EXTRA_DIST = \
	multisphinx.i \
	includes_java.i \
	includes_python.i \
	common.i 

libmultisphinx_jni.so: multisphinx_wrap_java.o
	$(CC) -shared -o $@ multisphinx_wrap_java.o \
		-L../src/libsphinxbase/.libs -L../src/libpocketsphinx/.libs \
		-lsphinxbase -lpocketsphinx

multisphinx_wrap_java.c: multisphinx.i common.i pocketsphinx_java.i
	mkdir -p edu/cmu/multisphinx
	swig -I../include -I../../sphinxbase/include \
		-java -package edu.cmu.multisphinx \
		-o $@ -outdir edu/cmu/multisphinx multisphinx.i

clean:
	$(RM) -r *.so *_wrap*.c *_wrap*.o *.so \
	multisphinx.py build *.egg-info *.class edu/cmu/*/*.class
