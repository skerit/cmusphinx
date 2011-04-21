#!/usr/bin/env python

"""
Simulated arc buffers for python prototyping.
"""

import numpy as np

class Arc(object):
    """
    Arc in a simulated arc buffer.
    """
    __slots__ = [ "word", "sf", "ef", "pscr", "lscr", "rc" ]

class ArcBuffer(list):
    """
    Simulated arc buffer class.
    """
    def __init__(self, arcfile):
        self.read_arcs(arcfile)

    def read_arcs(self, arcfh):
        if not isinstance(arcfh, file):
            arcfh = file(arcfh)
        cur_frame = -1
        for spam in arcfh:
            fields = spam.strip().split()
            if fields[0] == "__START_FRAME":
                cur_frame = int(fields[1])
            elif fields[0] == "__END_FRAME":
                cur_frame = -1
            else:
                arc = Arc()
                arc.word, arc.sf, arc.ef, arc.pscr, arc.lscr = fields[0:5]
                rc = tuple(tuple(map(int, x.split(":"))) for x in fields[5:])
                arc.rc = np.zeros(max(x[0] for x in rc) + 1)
                arc.rc.put(*zip(*rc))
                self.append(arc)
                
if __name__ == "__main__":
    import multisphinx
    m = multisphinx.Mdef(None, "../../data/hub4wsj_sc_8k/mdef")
    d = multisphinx.Dict(m, "../../data/bcb05cnp.dic")
    d2p = multisphinx.DictToPid(m, d)
    a = ArcBuffer("00000000.arc")
    for arc in a:
        wid = d.wordid(arc.word)
        rcmap = d2p.get_rcmap(wid)
        print arc.word, arc.sf, arc.ef,
        for i in range(m.n_ciphone()):
            if rcmap[i] < len(arc.rc) and not m.is_fillerphone(i):
                print "%s:%d" % (m.ciphone_str(i), arc.rc[rcmap[i]]),
        print
