'''
Created on 2011-04-07

@author: dhuggins
'''

import multisphinx

if __name__ == '__main__':
    sf = multisphinx.SearchFactory(["-lm", "../test/data/bn10000.3g.arpa",
                                "-hmm", "../test/data/hub4wsj_sc_8k",
                                "-dict", "../test/data/bn10000.dic", "-samprate", "11025"])
    sa = sf.create("state_align")
    al = multisphinx.Alignment(sf.d2p())
    al.add_words(["<s>", "A", "LONG", "DISCUSSION", "ABOUT", "HOW",
                  "WE", "WANT", "TO", "MAKE", "IT", "EASY", "FOR",
                  "PEOPLE", "TO", "BLEEP", "THINGS", "OUT", "</s>"])
    al.populate_ci()
    sa.set_alignment(al)
    sa.run()
    fb = sf.featbuf()
    foo = file("../test/data/chan3.raw")
    fb.producer_start_utt("chan3");
    fb.producer_process_raw(foo.read())
    fb.producer_end_utt()
    d = sf.dict()
    for w in al.words():
        print d.wordstr(w.id.wid), w.start, w.duration
    print sa.hyp()
    fb.producer_shutdown()