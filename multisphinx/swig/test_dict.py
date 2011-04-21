import multisphinx

m = multisphinx.Mdef(None, "../test/data/hub4wsj_sc_8k/mdef")
d = multisphinx.Dict(m, "../test/data/cmu07a.dic")
d2p = multisphinx.DictToPid(m, d)

wid = d.wordid("zwicky")
print wid
print d.wordstr(wid)
print m.ciphone_str(d.last_phone(wid))

xwdssid = d2p.rssid(d.last_phone(wid), d.second_last_phone(wid))
for i in range(0, m.n_ciphone()):
    print xwdssid.cimap[i],
print
for i in range(0, m.n_ciphone()):
    print xwdssid.ssid[i],
print
print xwdssid.n_ssid

