"""
purpose: convert native photon field data to a unified format.
These files have to be stored at share/crpropa/Scaling/
for CRPropa to read. The 10 supported photon fields are:

	IRB_Dominguez11, IRB_Finke10, IRB_Franceschini08, IRB_Gilmore12, IRB_Kneiske04,
	IRB_Stecker05, IRB_Stecker16_upper, IRB_Stecker16_lower, URB_Protheroe96, CMB

output: .txt files named by their field name. Format of output file is:
	-> lines starting with "#": ignored info lines
	-> line 1: photon energy e in eV
	-> line 2: redshift z at which the field is defined
	-> line 3+: field density d in 1/eVcm³: d(e0,z0), ... ,d(e0,zn)
												:	   %	  :
											d(en,z0), ... ,d(en,zn)
Written by Mario Hörbe (mario.hoerbe@rub.de / mario.hoerbe@physics.ox.ac.uk)
"""

import numpy as np

eV = 1.60217657e-19  # [J]
cm3 = 1e-6  # [m³]


def writeField(photonField, eps=[]):
	with open(photonField.name + ".txt", 'w') as f:
		f.write("# " + photonField.info + "\n")
		
		# photon energy in which field is defined
		f.write("# photon energy / eV:\n")
		if len(eps) == 0:
			eps = photonField.data[0][0]
		np.savetxt(f, [eps / eV], fmt="%1.6e", delimiter=" ")

		# redshift where field is defined
		f.write("# redshift:\n")
		energySlices = []
		if "redshift" not in locals():  # treat isotropy as large redshift range
			redshift = [0, 10]
			for z in redshift:
				energySlices.append(photonField.getDensity(eps) * eV * cm3)
		else:
			redshift = photonField.redshift
			for z in redshift:
				energySlices.append(photonField.getDensity(eps, z) * eV * cm3)
		np.savetxt(f, [redshift], fmt="%1.2f", delimiter=" ")
		
		# field's photon density
		f.write("# photon field density / 1/eVcm^3:\n")
		np.savetxt(f, list(map(list, zip(*energySlices))), fmt="%1.6e", delimiter=" ")
	print("done: " + photonField.name)


if __name__ == "__main__":

	import photonField

	writeField(photonField.EBL_Dominguez11())
	writeField(photonField.EBL_Finke10())
	writeField(photonField.EBL_Franceschini08())
	writeField(photonField.EBL_Gilmore12())
	writeField(photonField.EBL_Kneiske04())
	writeField(photonField.EBL_Stecker05())
	writeField(photonField.EBL_Stecker16("upper"))
	writeField(photonField.EBL_Stecker16("lower"))

	URB = photonField.URB_Protheroe96()
	eps = np.logspace(np.log10(URB.getEmin()), np.log10(URB.getEmax()), 100)
	writeField(URB, eps)

	CMB = photonField.CMB()
	eps = np.logspace(np.log10(CMB.getEmin()), np.log10(CMB.getEmax()), 100)
	writeField(CMB, eps)
