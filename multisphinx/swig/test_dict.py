import multisphinx

m = multisphinx.Mdef(None, "../test/data/hub4wsj_sc_8k/mdef")
d = multisphinx.Dict(m, "../test/data/cmu07a.dic")
d2p = multisphinx.DictToPid(m, d)

print d.wordid("zwicky")
print d.wordstr(d.wordid("zwicky"))
print m.ciphone_str(d.last_phone(d.wordid("zwicky")))

rcmap = d2p.get_rcmap(d.wordid("zwicky"))
for i in range(0, m.n_ciphone()):
    print rcmap[i],
print
