import numpy as np

with open('original_vmd.xyz', 'w') as f:
    with open('crys.crds', 'r') as g:
        crds = np.genfromtxt(g)
    f.write('{:}\n\n'.format(crds.shape[0]))
    no_si = int(crds.shape[0]/3)
    for i in range(no_si):
        f.write('Si{:}    {:.8f}      {:.8f}      {:.8f}\n'.format(i, crds[i,0], crds[i,1], crds[i,2]))
    for i in range(no_si, crds.shape[0]):
        f.write('O{:}    {:.8f}      {:.8f}      {:.8f}\n'.format(i-no_si, crds[i,0], crds[i,1], crds[i,2]))

    
