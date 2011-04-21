
typedef struct featbuf_s {
} FeatBuf;

%extend FeatBuf {
	FeatBuf(Config *c) {
		return featbuf_init(c);
	}
	~FeatBuf() {
		featbuf_free($self);
	}
	void producer_start_utt(char *uttid) {
		/* FIXME: Need to copy uttid or something. */
		/* FIXME: Also, exceptions, etc. */
		featbuf_producer_start_utt($self, uttid);
	}
	void producer_end_utt() {
		/* FIXME: exceptions, etc. */
		featbuf_producer_end_utt($self);
	}
	void producer_shutdown() {
		featbuf_producer_shutdown($self);
	}
	int producer_process_raw(short const *SDATA, size_t NSAMP, bool full_utt=false) {
		return featbuf_producer_process_raw($self, SDATA,
						    NSAMP / sizeof(int16),
						    full_utt);
	}
%pythoncode %{
def producer_process_audio(self, audio, full_utt=False):
    return self.producer_process_raw(audio.astype("int16").tostring(), full_utt)
%}
	int producer_process_cep(void *indata, int inlen, bool full_utt=false) {
		int ceplen = feat_cepsize(featbuf_get_fcb($self));
		return featbuf_producer_process_cep($self, indata,
						    inlen / (ceplen * sizeof(mfcc_t)),
						    full_utt);
	}
%pythoncode %{
def producer_process_mfcc(self, mfcc, full_utt=False):
    return self.producer_process_raw(mfcc.ravel().astype("float32").tostring(),
				     full_utt)
%}
};

