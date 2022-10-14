#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     17/01/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import ddatoolbox
#
path = "./ShimadzuLCMS9040mzxml/*"
#
# m/z value of target metabolite
#
#targetmz = 718.5381 # PE (34:1) [M+H]+ , C39H76NO8P.
targetmz = 805.4867 # PI 32:2 [M-H]-, C41H75O13P
#
# m/z threshold for dot product calclation & product ion clustering
#
threshold = 0.01
#
# Data loading. How supperting Shimadzu Q-TOF mzxml format.
#
msmsdata = ddatoolbox.loadxmlflles(path, targetmz)
#
# Construction of clieques
#
DDAcliques = ddatoolbox.constructonofDDAclique(msmsdata, threshold = 0.01, minimumsize = 3)
#
# Generation of averaged spectra for each clique
#
for i, msmslist in enumerate(DDAcliques):
    #
    # Header
    #
    print("###")
    print(str(i) + " th clique")
    print("total number of " + str(len(msmslist)) + " spectra included")
    #
    # Mean retention time and precursor m/z
    #
    mean_rt = ddatoolbox.calcrt(msmslist, msmsdata)
    mean_precursormz = ddatoolbox.calcmz(msmslist, msmsdata)
    print("Retention time: {:.1f}".format(mean_rt))
    print("Prec m/z: {:.5f}".format(mean_precursormz))
    #
    # Generation of averaged spectra
    #
    averaged_spectra_hash = ddatoolbox.constructonofProductIonclique(msmslist, msmsdata, threshold)
    for fragmentmz, mean_intensity in averaged_spectra_hash.items():
        print("\t {:.5f} {:>4d}".format(fragmentmz, mean_intensity))









