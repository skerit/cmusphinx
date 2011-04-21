CFLAGS := -g -Wall -DPIC -fPIC
CPPFLAGS := -I/usr/lib/jvm/default-java/include -I../include

all: libpocketsphinx_jni.so

libpocketsphinx_jni.so: pocketsphinx_wrap.o multisphinx_wrap.o
	$(CC) -shared -o $@ pocketsphinx_wrap.o \
		-L../src/libsphinxbase/.libs -L../src/libpocketsphinx/.libs \
		-lsphinxbase -lpocketsphinx

pocketsphinx_wrap.c: pocketsphinx.i
	mkdir -p edu/cmu/pocketsphinx
	swig -I../include -I../../sphinxbase/include \
		-java -package edu.cmu.pocketsphinx \
		-outdir edu/cmu/pocketsphinx pocketsphinx.i

multisphinx_wrap.c: multisphinx.i
	mkdir -p edu/cmu/multisphinx
	swig -I../include -I../../sphinxbase/include \
		-java -package edu.cmu.multisphinx \
		-outdir edu/cmu/multisphinx multisphinx.i

clean:
	$(RM) *.so *_wrap.c *_wrap.o *.class edu/cmu/*/*.class
