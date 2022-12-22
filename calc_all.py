import time
import subprocess
import photonField
import calc_elasticscattering as es
import calc_electromagnetic as em
import calc_pairproduction as bh
import calc_photodisintegration as pdi
import calc_photopionproduction as ppp
import calc_synchrotron as syn
import calc_photonFields as cpf

def nuclear_decay():
    """Creates nuclear decay tables"""

    print("#"*50)
    print("Calculate nuclear decay tables.\n")
    t1 = time.time()
    import calc_decay
    t2 = time.time()
    print("\nNuclear decay tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def elastic_scattering(fields: list):
    """Creates elastic scattering tables
    
    Input
        fields: list of CRPropa photonField class instances
    """

    print("#"*50)
    print("Calculate elastic scattering.\n")
    t1 = time.time()
    for field in fields:
        print(field.name)
        es.process(field)
    t2 = time.time()
    print("\nElastic scattering tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")


def EM_processes(fields: list):
    """Calculate electromagnetic processes
    
    Input
        fields: list of CRPropa photonField class instances
    """

    t1 = time.time()
    print("#"*50)
    print("Calculate electromagnetic processes:")
    print("\t- Pair production")
    print("\t- Double pair production")
    print("\t- Triple pair production")
    print("\t- Inverse Compton scattering")
    for field in fields:
        print(field.name)
        em.process(em.sigmaPP, field, 'EMPairProduction')
        em.process(em.sigmaDPP, field, 'EMDoublePairProduction')
        em.process(em.sigmaTPP, field, 'EMTripletPairProduction')
        em.process(em.sigmaICS, field, 'EMInverseComptonScattering')
    t2 = time.time()
    print("\Electromagnetic processes generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n") 

def nuclear_mass():
    """Creates nuclear mass tables"""

    print("#"*50)
    print("Calculate nuclear mass tables.\n")
    t1 = time.time()
    import calc_mass
    t2 = time.time()
    print("\nMass tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def BH_pair_production(fields: list):
    """Creates hadronic pair production
    
    Input
        fields: list of CRPropa photonField class instances
    """

    print("#"*50)
    print("Calculate Bethe-Heitler pair production.\n")
    t1 = time.time()
    bh.reformat_secondary_rates()
    for field in fields:
        print(field.name)
        bh.process(field)
    t2 = time.time()
    print("\nBH pair production tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def photo_disintegration(rateFields: list, emissionFields: list):
    """Calculate photo disintegration
    
    Input
        fields: list of CRPropa photonField class instances
    """

    print("#"*50)
    print("Calculate photo disintegration.\n")
    t1 = time.time()
    for field in rateFields:
        print(field.name)
        pdi.processRate(field)
    for field in emissionFields:
        print(field.name)
        pdi.processEmission(field)
    t2 = time.time()
    print("\nPhoto disintegration tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

#TODO That needs to be updated to work with arbitrary photon field classes.
def photon_fields():
    """Photon field processing
    
    This calls routines to calculate, e.g., the redshift scaling files.
    """

    print("#"*50)
    print("Process Photon fields\n")
    t1 = time.time()
    cpf.process()
    t2 = time.time()
    print("\nPhotonfield data generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def photopion_production(fields: list):
    """Creates photo pion production tables
    
    Input
        fields: list of CRPropa photonField class instances
    """

    print("#"*50)
    print("Calculate photo pion production\n")
    t1 = time.time()
    for field in fields:
        print(field.name)
        ppp.process(field)
    t2 = time.time()
    print("\nPhoto pion production tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def synchrotron():
    """Creates synchrotron spectra"""

    print("#"*50)
    print("Calculate synchrotron radiation\n")
    t1 = time.time()
    syn.process()
    t2 = time.time()
    print("\nSynchrotron tables generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def compress():
    """Compressing the data directory into a tar.gz-ball
    
    The compressed files will be named 'data-YYYY-MM-DD.tar.gz'
    """

    print("#"*50)
    print("Compressing the ./data directory.")
    t1 = time.time()
    datestr = "-".join([str(x) for x in time.localtime()[:3]])
    subprocess.run(["tar", "-czf", "data-"+datestr+".tar.gz", "./data"])
    t2 = time.time()
    print("\nCompressed files generated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def calc_checksum():
    """Calculate the md5-checksum of the tarball
    
    The Checksum is stored in a file called 
    data-YYYY-MM-DD.tar.gz-CHECKSUM
    """

    print("#"*50)
    print("Calculating the checksum.")
    t1 = time.time()
    datestr = "-".join([str(x) for x in time.localtime()[:3]])
    checksum = subprocess.run(["md5sum", "data-"+datestr+".tar.gz"],
                    capture_output=True, text=True).stdout
    with open("data-"+datestr+".tar.gz-CHECKSUM", 'w') as f:
        f.write(checksum)
    t2 = time.time()
    print("\nChecksum calculated in {} seconds.".format(round(t2-t1, 2)))
    print("#"*50+"\n")

def createPhotonTargetInteractions(fields: list):
    """Create all interaction tables that depend on a photon field
    
    Input
        fields: list of list of CRPropa photonField class instances
    """
    elastic_scattering(fields)
    EM_processes(fields)
    BH_pair_production(fields)
    photo_disintegration(fields, fields)
    photopion_production(fields)

    #TODO: Include the creation of the scaling files

def createCRPropaDefault():
    """Creating and compressing all default files 
    needed for the default CRPropa installation.

    Note: The default CRPropa data set uses for some
    interactions (elastic scattering, Bethe-Heitler pair production)
    and the emission files of photo-disintegration only
    a subset of all default photon fields. This is done
    to reduce the amount of data that need to be shipped 
    with the code.
    """

    reduced_fields = [
        photonField.CMB(),
        photonField.EBL_Gilmore12(),
        photonField.URB_Protheroe96()
        ]

    fields_cmbebl = [
        photonField.CMB(),
        photonField.EBL_Kneiske04(),
        photonField.EBL_Stecker05(),
        photonField.EBL_Franceschini08(),
        photonField.EBL_Finke10(),
        photonField.EBL_Dominguez11(),
        photonField.EBL_Gilmore12(),
        photonField.EBL_Stecker16('lower'),
        photonField.EBL_Stecker16('upper')
    ]

    fields_urb = [
        photonField.URB_Protheroe96(),
        photonField.URB_Fixsen11(),
        photonField.URB_Nitu21()
    ]

    nuclear_decay()
    nuclear_mass()
    elastic_scattering(reduced_fields)
    EM_processes(fields_cmbebl+fields_urb)
    BH_pair_production(fields_cmbebl)
    photo_disintegration(fields_cmbebl+fields_urb, reduced_fields)
    photon_fields()
    photopion_production(fields_cmbebl+fields_urb)
    synchrotron()
    compress()
    calc_checksum() 


if __name__ == "__main__":
    createCRPropaDefault()

    
    
    
    
    
    
    
    
    
    
    
