Mikko = open('g_eff.itx')
Mikkodat = Mikko.read().split('\n')
Mikkodat = Mikkodat[7:]
ge_eff = []
gs_eff = []
T = []
for i,v in enumerate(Mikkodat):
    T.append(float(v.split()[0])*10**6) #in eV
    ge_eff.append(float(v.split()[7]))
    gs_eff.append(float(v.split()[8]))
ge_eff = np.array(ge_eff)
gs_eff = np.array(gs_eff)
T = np.array(T)
@function
