#-------------------------------------------------------------------------------
# Name:        ddatoolbox
# Purpose:
#
# Author:      FumioMatsuda
#
# Created:     07/10/2022
# Copyright:   (c) FumioMatsuda 2022
# Licence:     MIT lisence
#-------------------------------------------------------------------------------

import xml.etree.ElementTree as ET
import base64, struct,itertools, glob, re
import numpy as np


def dot_product(mshash1, mshash2, threshold = 0.01):
    """
    A function for calclation of dot product score.
    mshash1 : Spec hash 1
    mshash2 : Spec hash 2
    threshold: Threshold for m/z comparison
    """
    norm1 = mshash1['norm']
    norm2 = mshash2['norm']
    prec1 = mshash1['precursor_mz']
    prec2 = mshash2['precursor_mz']
    if abs(prec1 - prec2) > threshold:
        return(0.0)

    dotproduct = 0.0
    for (precmz1, precmz2) in itertools.product(mshash1['spectrum'].keys(), mshash2['spectrum'].keys()):
        if abs(float(precmz1) - float(precmz2)) < threshold:
            dotproduct = dotproduct + (mshash1['spectrum'][precmz1] / norm1) * (mshash2['spectrum'][precmz2] / norm2)
    return(dotproduct)

def decodeBase64(base64text):
    """
    A function to decode Base64 derived from
    https://groups.google.com/g/spctools-discuss/c/qK_QThoEzeQ
    """
    decoded = base64.b64decode(base64text)
    tmp_size = len(decoded)/8
    unpack_format1 = ">%dQ" % tmp_size

    idx = 0
    mz_list = []
    intensity_list = []

    for tmp in struct.unpack(unpack_format1,decoded):
        #print(tmp)
        tmp_i = struct.pack("Q",tmp)
        tmp_f = struct.unpack("d",tmp_i)[0]
        if( idx % 2 == 0 ):
            mz_list.append( float(tmp_f) )
        else:
            intensity_list.append( float(tmp_f) )
        idx += 1
    basepeakintensity = max(intensity_list)
    mz_list_selected = []
    intensity_list_selected = []
    for mz, intensity in zip(mz_list, intensity_list):
        if intensity/ basepeakintensity * 1000 < 5:
            continue
        mz_list_selected.append(mz)
        intensity_list_selected.append(intensity)


    norm = np.linalg.norm(intensity_list_selected, ord=2)
    return(norm, dict(zip(mz_list_selected, intensity_list_selected)))


def constructonofDDAclique(msmsdata, threshold = 0.01, minimumsize = 3):
    """
    Construction of DDA cliques
    This version does not use xml.etree.ElementTree
    """
    mzmzlist = msmsdata.keys()
    cliques = []
    for msms1 in sorted(mzmzlist, key=lambda x: msmsdata[x]['basePeakIntensity'], reverse = True):
        hitclique = []
        hitclique_averagescore = []
        for clique in cliques:
            averagescore_temp = []
            for msms2 in clique:
                dotproduct = dot_product(msmsdata[msms1], msmsdata[msms2], threshold = threshold)
                averagescore_temp.append(dotproduct)
                if not dotproduct > 0.90:
                    break
            else:
                hitclique.append(clique)
                hitclique_averagescore.append(np.mean(averagescore_temp))
        if hitclique == []:
            cliques.append([msms1])
            continue

        hitclique_order_sorted = sorted(range(len(hitclique_averagescore)), key=lambda x: len(hitclique[x]) + hitclique_averagescore[x], reverse = True)
        clique = hitclique[hitclique_order_sorted[0]]
        clique.append(msms1)
        #print(msms1, clique)
    cliquelist_sorted = sorted(cliques, key=lambda x: len(x), reverse = True)
    mzmzlistset = set(mzmzlist)

    cliquelist_selected = []
    for clique in cliquelist_sorted :
        if len(clique) < minimumsize:
            continue
        cliquelist_selected.append(clique)
    return(cliquelist_selected)


def constructonofProductIonclique(msmslist, msmsdata, threshold):
    """
    Construction of product ion cliques
    This version does not use xml.etree.ElementTree

    """
    cliques = []
    for msms1 in msmslist:
        msms1_precs = list(msmsdata[msms1]['spectrum'].keys())
        for msms1_prec in msms1_precs:
            label1 = msms1+":"+str(msms1_prec)
            #print(msms1)
            hitclique = []
            hitclique_averagescore = []
            for clique in cliques:
                averagescore_temp = []
                for label2 in clique:
                    (msms2, msms2_prec) = label2.split(":")
                    diff = abs(msms1_prec - float(msms2_prec))
                    averagescore_temp.append(diff)
                    if not diff < threshold:
                        break
                else:
                    hitclique.append(clique)
                    hitclique_averagescore.append(np.mean(averagescore_temp))
            if hitclique == []:
                cliques.append([label1])
                continue

            hitclique_order_sorted = sorted(range(len(hitclique_averagescore)), key=lambda x: len(hitclique[x]) - hitclique_averagescore[x], reverse = True)

            clique = hitclique[hitclique_order_sorted[0]]
            clique.append(label1)

    msmslist_number = len(msmslist)
    temp_hash = {}
    for msmsclique in cliques:
        msmsclique_number = len(msmsclique)

        if msmsclique_number/msmslist_number > 0.7:
            fragmentmzlist = []
            fragmentintensitylist = []
            precmzlist = []
            for fragment in msmsclique:
                (nodeid, fragmentmz) = fragment.split(":")
                #print(nodeid, fragmentmz, msmsdata[nodeid]['spectrum'].keys())
                fragmentmzlist.append(float(fragmentmz))
                intensity = msmsdata[nodeid]['spectrum'][float(fragmentmz)]
                basePeakIntensity = float(msmsdata[nodeid]['basePeakIntensity'])
                fragmentintensitylist.append(float(intensity)/basePeakIntensity)
            mean_mz = np.mean(fragmentmzlist)

            if len(fragmentmzlist) > 30:
                mean_mz = np.median(fragmentmzlist)


            mean_intensity = int(sum(fragmentintensitylist)/len(fragmentintensitylist) * 999)
            temp_hash[mean_mz] = mean_intensity
    return(temp_hash)



def loadxmlflles(path, targetmz, threshold = 0.02):
    """
    A function to load Shimadzu mzxml files at "path"
    targetmz: MS2 data obtained from precursor ions within targetmz+- threshold were collected.

    """
    msmsdata = {}
    files = glob.glob(path)


    for file in files:

        data_number = file.split('\\')
        data_number = data_number[-1]
        data_number = data_number.replace(".mzXML", "")

        print("loading files", file, data_number)


        tree = ET.parse(file)
        root = tree.getroot()
        for scan in root.iter('{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}scan'):

            if scan.attrib["msLevel"] != "2":
                continue
            spectrum_number = scan.attrib["num"]
            basePeakIntensity = scan.attrib['basePeakIntensity']
            retentionTime = scan.attrib['retentionTime']

            for peaks in scan.findall('{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}precursorMz'):
                precursormz = peaks.text

                lowerboundary = targetmz - threshold
                upperboundary = targetmz + threshold

                if lowerboundary < float(precursormz) < upperboundary:
                    break
            else:
                continue


            for peaks in scan.findall('{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}peaks'):
                if peaks.text == None:
                    continue
                norm, msmsdict = decodeBase64(peaks.text)
                #
                # Discurded if number of production ion is less than one
                #
                if len(msmsdict.keys()) >= 2:
                    msmsspecid = data_number+"_"+spectrum_number
                    msmsdata[msmsspecid] = {'precursor_mz': float(precursormz), 'norm':norm,'spectrum': msmsdict, 'retentionTime': retentionTime, 'basePeakIntensity':basePeakIntensity}
    return(msmsdata)

def calcrt(msmslist, msmsdata):
    """
    A function to calculate mean retention time of given msmslist in msmsdata
    """
    rtlist = []
    for msmsid in sorted(msmslist):
        temprt =msmsdata[msmsid]['retentionTime']
        temprt = re.sub("PT", "", temprt)
        temprt = re.sub("S", "", temprt)
        rtlist.append(float(temprt))


    return(sum(rtlist)/len(rtlist))

def calcmz(msmslist, msmsdata):
    """
    A function to calculate mean precursor m/z of given msmslist in msmsdata
    """
    precursormzlist = []
    for msmsid in msmslist:
        precursormzlist.append(float(msmsdata[msmsid]['precursor_mz']))
    mean_precursormz = np.mean(precursormzlist)
    if len(precursormzlist) > 30:
        mean_precursormz = np.median(precursormzlist)

    return(mean_precursormz)

