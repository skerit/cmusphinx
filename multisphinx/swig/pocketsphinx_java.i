%include <arrays_java.i>

/* Special typemap for arrays of audio. */
%typemap(in) (short const *SDATA, size_t NSAMP) {
	$1 = (short const *) JCALL2(GetShortArrayElements, jenv, $input, NULL);
	$2 = JCALL1(GetArrayLength, jenv, $input);
};
%typemap(freearg) (short const *SDATA, size_t NSAMP) {
	JCALL3(ReleaseShortArrayElements, jenv, $input, $1, 0);
};
%typemap(jni) (short const *SDATA, size_t NSAMP) "jshortArray"
%typemap(jtype) (short const *SDATA, size_t NSAMP) "short[]"
%typemap(jstype) (short const *SDATA, size_t NSAMP) "short[]"
%typemap(javain) (short const *SDATA, size_t NSAMP) "$javainput"

