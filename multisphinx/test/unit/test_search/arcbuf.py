#!/usr/bin/env python

"""
Simulated arc buffers for python prototyping.
"""

import numpy as np
import multisphinx
import lattice

class Arc(object):
    """
    Arc in a simulated arc buffer.
    """
    __slots__ = [ "word", "sf", "ef", "pscr", "lscr", "rcdelta" ]

class ArcBuffer(list):
    """
    Simulated arc buffer class.
    """
    def __init__(self, arcfile, d2p):
        self.d2p = d2p
        self.read_arcs(arcfile)

    def read_arcs(self, arcfh):
        if not isinstance(arcfh, file):
            arcfh = file(arcfh)
        cur_frame = -1
        self.n_frames = -1
        n_ci = self.d2p.get_mdef().n_ciphone()
        dic = self.d2p.get_dict()
        for spam in arcfh:
            fields = spam.strip().split()
            if fields[0] == "__START_FRAME":
                cur_frame = int(fields[1])
            elif fields[0] == "__END_FRAME":
                cur_frame = -1
            else:
                arc = Arc()
                arc.word, arc.sf, arc.ef, arc.pscr, arc.lscr = fields[0:5]
                arc.sf = int(arc.sf)
                arc.ef = int(arc.ef) + 1 # Adjust to be non-inclusive
                self.n_frames = max(self.n_frames, arc.ef)
                arc.pscr = int(arc.pscr)
                arc.lscr = int(arc.lscr)
                rc = dict(map(int, x.split(":")) for x in fields[5:])
                arc.rcdelta = np.ones(n_ci, 'uint16')
                rcmap = d2p.get_rcmap(dic.wordid(arc.word))
                for i in xrange(n_ci):
                    arc.rcdelta[i] = rc.get(rcmap[i], 65535)
                self.append(arc)
        self.eflist = self[:]
        self.eflist.sort(key=lambda x:x.ef)
        self.efidx = np.ones(self.n_frames + 1, 'i') * -1
        cur_frame = -1
        for i, arc in enumerate(self.eflist):
            if arc.ef > cur_frame:
                self.efidx[arc.ef] = i
                cur_frame = arc.ef

    def ending_at(self, ef):
        start = self.efidx[ef]
        if start == -1:
            return
        for arc in self.eflist[start:]:
            if arc.ef != ef:
                return
            yield arc

def build_arcbuf_lattice(arcbuf):
    dag = lattice.Dag()
    dag.nodelist = []
    dic = arcbuf.d2p.get_dict()
    def build_arc_node(word, frame):
        node = dag.Node(word, frame)
        dag.add_node(node)
        if frame == 0:
            node.score = 0
        else:
            for incoming in arcbuf.ending_at(node.frame):
                src = dag.find_node(incoming.word, incoming.sf)
                rcdelta = incoming.rcdelta[dic.first_phone(dic.wordid(incoming.word))]
                if rcdelta == 65535:
                    continue
                incoming_score = incoming.pscr - rcdelta
                link = dag.Link(src, node, incoming_score
                                - incoming.lscr - src.score)
                src.exits.append(link)
                node.entries.append(link)
                node.score = max(node.score, incoming_score)
        return node
    for outgoing in arcbuf:
        if outging.word == "</s>" and outgoing.ef != arcbuf.n_frames: continue
        if dag.find_node(outgoing.word, outgoing.sf): continue
        build_arc_node(outgoing.word, outgoing.sf)
    end_pscr = -999999999
    for incoming in arcbuf.ending_at(arcbuf.n_frames):
        if incoming.word != "</s>": continue
        if incoming.pscr > end_pscr:
            src = dag.find_node(incoming.word, incoming.sf)
            dag.end = src
            dag.final_ascr = incoming.pscr - incoming.lscr - src.score
    dag.start = dag.find_node("<s>", 0)
    dag.remove_unreachable()
    return dag
                
if __name__ == "__main__":
    m = multisphinx.Mdef(None, "../../data/hub4wsj_sc_8k/mdef")
    d = multisphinx.Dict(m, "../../data/bcb05cnp.dic")
    d2p = multisphinx.DictToPid(m, d)
    a = ArcBuffer("00000000.arc", d2p)
    dag = build_arcbuf_lattice(a)
    dag.dag2dot("00000000.dot", lambda x:x.ascr)
