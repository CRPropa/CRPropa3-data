import os
import numpy as np
from scipy.integrate import trapz
import photonField
import gitHelp as gh

# Calculate the integral I(z) of the EBL spectral number density dn/depsilon as function of redshift z.
# The ratio s(z) = I(z)/I(0) serves as a global scaling factor for all interactions with the EBL.
# In contrast to CRPropa 2, the photon spectrum is integrated over the whole tabulated range.

# Note: Redshift Scaling is now natively supported in CRPropa's PhotonBackground class. Therefore, 
# the scaling files are not included in the download data and this script is only kept for comparison.

fields = [
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12(),
    photonField.EBL_Stecker16('upper'),
    photonField.EBL_Stecker16('lower')]

for field in fields:
    print(field.name)

    data = field.data  # dictionary {z : (eps, n(eps))}
    tz = np.array(list(data.keys()))
    tz.sort()
    ts = np.zeros_like(tz)

    for i, z in enumerate(tz):
        eps, n = data[z]
        ts[i] = trapz(n, eps)

    ts /= ts[0]
    
    # output folder
    folder = 'data/Scaling'
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    try:
        git_hash = gh.get_git_revision_hash()
        header = 'Produced with crpropa-data version: '+git_hash+'\nredshift\t global evolution factor'
    except:
        header = 'redshift\t global evolution factor'
    np.savetxt(
        'data/Scaling/%s_scaling.txt' % field.name,
        np.c_[tz, ts],
        fmt='%.2f\t%.4e',
        header=header)
