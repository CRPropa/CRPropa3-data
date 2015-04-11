import numpy as np
from scipy.integrate import trapz
import photonField


# Calculate the integral I(z) of the EBL spectral number density dn/depsilon as function of redshift z.
# The ratio s(z) = I(z)/I(0) serves as a global scaling factor for all interactions with the EBL.
# Contrary to CRPropa 2, the photon spectrum is integrated over the whole tabulated range.
def process(ebl, fname):
    data = ebl.data  # dictionary {z : (eps, n(eps))}

    tz = np.array(data.keys())
    tz.sort()
    ts = np.zeros_like(tz)

    for i, z in enumerate(tz):
        eps, n = data[z]
        ts[i] = trapz(n, eps)

    ts /= ts[0]
    np.savetxt(fname, np.c_[tz, ts], fmt='%.2f\t%.4e',
        header='redshift\t global evolution factor')


process(photonField.EBL_Kneiske04(), 'data/scaling_Kneiske04.txt')
process(photonField.EBL_Stecker05(), 'data/scaling_Stecker05.txt')
process(photonField.EBL_Finke10(),   'data/scaling_Finke10.txt')
process(photonField.EBL_Franceschini08(), 'data/scaling_Franceschini08.txt')
process(photonField.EBL_Dominguez11(), 'data/scaling_Dominguez11.txt')
process(photonField.EBL_Gilmore12(), 'data/scaling_Gilmore12.txt')

