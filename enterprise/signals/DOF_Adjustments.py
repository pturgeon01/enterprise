from enterprise.signals.parameter import function
import numpy as np
import urllib.request
url = "https://raw.githubusercontent.com/pturgeon01/enterprise/master/enterprise/signals/g_eff.itx"
urllib.request.urlretrieve(url, "g_eff.itx")

File = open('g_eff.itx')
Filedat = File.read().split('\n')
Filedat = Filedat[7:]
ge_eff = []
gs_eff = []
T = []
for i,v in enumerate(Filedat):
    d = v.split()
    T.append(float(d[0])*10**6) #in eV
    ge_eff.append(float(d[7]))
    gs_eff.append(float(d[8]))
T = np.array(T)
ge_eff = np.array(ge_eff)
gs_eff = np.array(gs_eff)

@function
def return_DOFge(input_value):
    # Find the index of the closest value in T to input_value
    input_value = np.array(input_value, ndmin=1)
    closest_indices = np.zeros(len(input_value), dtype=int)
    for i,v in enumerate(input_value):
        if np.abs(T[i] - v) < 1e40:
            closest_indices[i] = np.argmin(np.abs(T - v))
        else:
            closest_indices[i] = 1710 
    return ge_eff[closest_indices]

@function
def return_DOFgs(input_value):
    # Find the index of the closest value in T to input_value
    input_value = np.array(input_value, ndmin=1)
    closest_indices = np.zeros(len(input_value), dtype=int)
    for i,v in enumerate(input_value):
        if np.abs(T[i] - v) < 1e40:
            closest_indices[i] = np.argmin(np.abs(T - v))
        else:
            closest_indices[i] = 1710 
    return ge_eff[closest_indices]

