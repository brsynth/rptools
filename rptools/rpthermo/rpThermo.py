# import libsbml
# import argparse
# import sys #exit using sys exit if any error is encountered
# import os

# import io
# #import zipfile
# import tarfile
# import glob
# import tempfile
# import shutil

import logging
from equilibrator_api                import ComponentContribution, Q_
from rptools.rplibs                  import rpSBML
from rptools.rpthermo.rpEquilibrator import pathway


def runThermo(inFile, outFile,
              pathway_id='rp_pathway',
              ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15,
              logger=logging.getLogger(__name__)):
    """Given a tar input file, perform thermodynamics analysis for each rpSBML file.

    :param inFile: The path to the input file
    :param outFile: The path to the output file
    :param pathway_id: The id of the heterologous pathway of interest
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
    :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

    :type inFile: str
    :type outFile: str
    :type pathway_id: str
    :type tmpOutputFolder: str
    :type ph: float
    :type ionic_strength: float
    :type pMg: float
    :type temp_k: float
    :type stdev_factor: float

    :rtype: bool
    :return: Success or failure of the function
    """


    cc = initThermo(ph, ionic_strength, pMg, temp_k)

    rpsbml = rpSBML(inFile)

    results = pathway(rpsbml=rpsbml, cc=cc, pathway_id=pathway_id, update_rpsbml=True, logger=logger) # ignore the results since written to SBML file
 
    rpsbml.writeSBML(outFile)

    # #mnx_default_conc = json.load(open('data/mnx_default_conc.json', 'r'))
    # mnx_default_conc = json.load(open(os.path.join(os.path.dirname(os.path.abspath( __file__ )), 'data', 'mnx_default_conc.json'), 'r'))

    return True


# used to initialise and download the data for equilibrator
def initThermo(ph, ionic_strength, pMg, temp_k):
    cc = ComponentContribution()
    cc.p_h = Q_(ph)
    cc.ionic_strength = Q_(str(ionic_strength)+' mM')
    cc.p_mg = Q_(pMg)
    cc.temperature = Q_(str(temp_k)+' K')

    return cc

# def runMDF_hdd(inputTar, outputTar, pathway_id='rp_pathway', thermo_id='dfG_prime_o', fba_id='fba_obj_fraction', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
#     """Given a tar input file, perform MDF analysis for each rpSBML file.

#     :param inputTar: The path to the input TAR file
#     :param outputTar: The path to the output TAR file
#     :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
#     :param thermo_id: The id of the thermodynamics id (Default: dfG_prime_o)
#     :param fba_id: The id of the FBA value (Default: fba_obj_fraction)
#     :param ph: The pH of the host organism (Default: 7.0)
#     :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
#     :param pMg: The pMg of the host organism (Default: 10.0)
#     :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
#     :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

#     :type inputTar: str
#     :type outputTar: str
#     :type pathway_id: str
#     :type tmpOutputFolder: str
#     :type ph: float
#     :type ionic_strength: float
#     :type pMg: float
#     :type temp_k: float
#     :type stdev_factor: float

#     :rtype: bool
#     :return: Success or failure of the function
#     """
#     rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)
#     with tempfile.TemporaryDirectory() as tmpInputFolder:
#         with tempfile.TemporaryDirectory() as tmpOutputFolder:
#             tar = tarfile.open(inputTar, mode='r')
#             tar.extractall(path=tmpInputFolder)
#             tar.close()
#             if len(glob.glob(tmpInputFolder+'/*'))==0:
#                 logging.error('Input file is empty')
#                 return False
#             for sbml_path in glob.glob(tmpInputFolder+'/*'):
#                 logging.debug('=========== '+str(sbml_path)+' ============')
#                 fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '') 
#                 rpsbml = rpSBML.rpSBML(fileName, path=sbml_path)
#                 rpequilibrator.rpsbml = rpsbml
#                 res = rpequilibrator.MDF(pathway_id, thermo_id, fba_id, stdev_factor, True) #ignore the results since written to SBML file
#                 #ignore res since we are passing write to SBML
#                 rpsbml.writeSBML(tmpOutputFolder)
#                 rpsbml = None
#             if len(glob.glob(tmpOutputFolder+'/*'))==0:
#                 logging.error('rpThermo has not produced any results')
#                 return False
#             with tarfile.open(outputTar, mode='w:gz') as ot:
#                 for sbml_path in glob.glob(tmpOutputFolder+'/*'):
#                     fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', ''))
#                     fileName += '.sbml.xml'
#                     info = tarfile.TarInfo(fileName)
#                     info.size = os.path.getsize(sbml_path)
#                     ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
#     return True


# def runEqSBtab_hdd(inputTar, outputTar, pathway_id='rp_pathway', fba_id=None, thermo_id='dfG_prime_o', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
#     """Given a tar input file, perform MDF analysis for each rpSBML file.

#     :param inputTar: The path to the input TAR file
#     :param outputTar: The path to the output TAR file
#     :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
#     :param fba_id: The id of the FBA value. Default sets all FBA values to 1.0 and if specified (Default: None)
#     :param thermo_id: The id of the thermodynamics id (Default: dfG_prime_o)
#     :param ph: The pH of the host organism (Default: 7.0)
#     :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
#     :param pMg: The pMg of the host organism (Default: 10.0)
#     :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
#     :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

#     :type inputTar: str
#     :type outputTar: str
#     :type pathway_id: str
#     :type tmpOutputFolder: str
#     :type ph: float
#     :type ionic_strength: float
#     :type pMg: float
#     :type temp_k: float
#     :type stdev_factor: float

#     :rtype: bool
#     :return: Success or failure of the function
#     """
#     rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)
#     with tempfile.TemporaryDirectory() as tmpInputFolder:
#         with tempfile.TemporaryDirectory() as tmpOutputFolder:
#             tar = tarfile.open(inputTar, mode='r')
#             tar.extractall(path=tmpInputFolder)
#             tar.close()
#             if len(glob.glob(tmpInputFolder+'/*'))==0:
#                 logging.error('Input file is empty')
#                 return False
#             for sbml_path in glob.glob(tmpInputFolder+'/*'):
#                 logging.debug('=========== '+str(sbml_path)+' ============')
#                 fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '') 
#                 rpsbml = rpSBML.rpSBML(fileName, path=sbml_path)
#                 rpequilibrator.rpsbml = rpsbml
#                 status = rpequilibrator.toNetworkSBtab(os.path.join(tmpOutputFolder, fileName+'.tsv'), pathway_id, thermo_id, fba_id, stdev_factor)
#                 rpsbml = None
#             if len(glob.glob(tmpOutputFolder+'/*'))==0:
#                 logging.error('rpThermo has not produced any results')
#                 return False
#             with tarfile.open(outputTar, mode='w:gz') as ot:
#                 for sbml_path in glob.glob(tmpOutputFolder+'/*'):
#                     fileName = str(sbml_path.split('/')[-1])
#                     info = tarfile.TarInfo(fileName)
#                     info.size = os.path.getsize(sbml_path)
#                     ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
#     return True
