CRPropa3-data
=============

Tools to generate data files for [CRPropa 3](https://github.com/CRPropa/CRPropa3).
An example to create the necessary data files for a custom photon field can be found in the [CRPropa3 documentation](https://crpropa.github.io/CRPropa3/index.html).


CRPropa's default data files
 - calc_all.py

Interactions between cosmic ray nuclei and background photons
 - calc_pairproduction.py : electron pair production
 - calc_photopionproduction.py : photo-pion production
 - calc_photodisintegration.py : photodisintegration
 - calc_elasticscattering.py   : elastic scattering on the nuclear structure

Interactions between cosmic ray photons / electrons and background photons
 - calc_electromagnetic.py
    - photon    : pair and double-pair production
    - electrons : triplet pair production and inverse Compton scattering

Other processes
 - calc_decay.py : nuclear decays
 - calc_synchrotron.py : synchrotron radiation of charged particles
 - calc_mass.py : table of nuclear masses

Helper modules
 - photonField.py     : collection of background photon fields (CMB, EBL, URB)
 - interactionRate.py : functions to calculate interaction rates with isotropic photon fields
 - gitHelp.py : function to include git-hash in tabulated data

Galactic magnetic lenses
 - create_lens.py : create a magnetic lens from backtracking simualtions. See
	 ./create_lens.py for details.
