#!/usr/bin/env python

import multisphinx

sf = multisphinx.SearchFactory(["-lm", "../test/data/bn10000.3g.arpa",
                                "-hmm", "../test/data/hub4wsj_sc_8k",
                                "-dict", "../test/data/bn10000.dic", "-samprate", "11025"])
ft = sf.create("fwdtree")
ff = sf.create(ft, "fwdflat")
ft.link(ff, "fwdtree")

fb = sf.featbuf()
foo = file("../test/data/chan3.raw")
ft.run()
ff.run()
fb.producer_start_utt("chan3")
while True:
    bar = foo.read(2048)
    if len(bar) == 0:
        break
    fb.producer_process_raw(bar)
fb.producer_end_utt()
print ft.hyp()
print ff.hyp()
for foo in ff.seg_iter():
    print foo
fb.producer_shutdown()