CFLAGS := -g -Wall -DPIC -fPIC
CPPFLAGS := -I/usr/lib/jvm/default-java/include -I../include $(shell python-config --cflags)

all: libmultisphinx_jni.so

libmultisphinx_jni.so: pocketsphinx_wrap_java.o multisphinx_wrap_java.o
	$(CC) -shared -o $@ pocketsphinx_wrap_java.o multisphinx_wrap_java.o \
		-L../src/libsphinxbase/.libs -L../src/libpocketsphinx/.libs \
		-lsphinxbase -lpocketsphinx

pocketsphinx_wrap_java.c: pocketsphinx.i common.i pocketsphinx_java.i
	mkdir -p edu/cmu/pocketsphinx
	swig -I../include -I../../sphinxbase/include \
		-java -package edu.cmu.pocketsphinx \
		-o $@ -outdir edu/cmu/pocketsphinx pocketsphinx.i

multisphinx_wrap_java.c: multisphinx.i common.i pocketsphinx_java.i
	mkdir -p edu/cmu/multisphinx
	swig -I../include -I../../sphinxbase/include \
		-java -package edu.cmu.multisphinx \
		-o $@ -outdir edu/cmu/multisphinx multisphinx.i

clean:
	$(RM) *.so *_wrap*.c *_wrap*.o *.class edu/cmu/*/*.class
