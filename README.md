CRPropa3-data
=============

Tools to generate data files for [CRPropa 3](https://github.com/CRPropa/CRPropa3).

Interactions between cosmic ray nuclei and background photons
 - calc_epp.py : electron pair production
 - calc_ppp.py : photo-pion production
 - calc_photodisintegration.py : photodisintegration
 - calc_elasticscattering.py   : elastic scattering on the nuclear structure
 An example to create the necessary data files for a custom photon field can be found in the [CRPropa3 documentation](https://crpropa.github.io/CRPropa3/index.html).

Interactions between cosmic ray photons / electrons and background photons
 - calc_electromagnetic.py
    - photon    : pair and double-pair production
    - electrons : triplet pair production and inverse Compton scattering

Other processes
 - calc_decay.py : nuclear decays
 - calc_synchrotron.py : synchrotron radiation of charged particles
 - calc_scaling : global redshift scaling of cosmic photon fields
 - calc_mass : table of nuclear masses

Helper modules
 - photonField.py     : collection of background photon fields (CMB, EBL, URB)
 - interactionRate.py : functions to calculate interaction rates with isotropic photon fields

Galactic magnetic lenses
 - create_lens.py : create a magnetic lens from backtracking simualtions. See
	 ./create_lens.py for details.
