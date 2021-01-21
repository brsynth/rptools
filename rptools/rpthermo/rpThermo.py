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

from rpEquilibrator import pathway

###################### Multi ###############

def chunkIt(seq, num):
    """Seperate an array into equal lengths

    :param seq: The arrray to seperate
    :param num: The number of chunks to seperate the array into

    :type seq: list
    :type num: int

    :rtype: list
    :return: 2D list of chunks
    """
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out


def singleThermo(sbml_paths, pathway_id, tmpOutputFolder, ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
    """Given a list of rpSBML input files, perform thermodynamics analysis. Less memory effecient than the _hdd method but faster

    :param sbml_paths: The list of rpSBML paths passed to calculate thermodynamics
    :param pathway_id: The id of the heterologous pathway of interest
    :param tmpOutputFolder: The path to the output folder to write the result rpSBML file
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
    :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

    :type sbml_paths: str
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
    rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k, stdev_factor=stdev_factor)
    for sbml_path in sbml_paths:
        logging.debug('Calculating the thermodynamics of the pathway '+str(pathway_id)+' for the file: '+str(sbml_path))
        file_name = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '')
        rpsbml = rpSBML.rpSBML(file_name, path=sbml_path)
        rpequilibrator.rpsbml = rpsbml
        res = rpequilibrator.pathway(pathway_id, True) #ignore the results because written in SBML file
        rpsbml.writeSBML(tmpOutputFolder)
        rpsbml = None
    return True

''' Seems like subprocess does not play well with multiprocessing
### concurent

import concurrent.futures

## Multiprocessing implementation of the thermodynamics package
#
#
def runThermo_multi_concurrent(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
    rpcache = rpCache.rpCache()
    cc_preprocess = rpcache.getCCpreprocess()
    kegg_dG = rpcache.getKEGGdG()
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
                jobs = {}
                #split the files "equally" between all workers
                for s_l in chunkIt(glob.glob(tmpInputFolder+'/*'), num_workers):
                    jobs[executor.submit(singleThermo, s_l, pathway_id, tmpOutputFolder, kegg_dG, cc_preprocess, ph, ionic_strength, pMg, temp_k)] = s_l
                for future in concurrent.futures.as_completed(jobs):
                    f_n = jobs[future]
                    try:
                        data = future.result()
                    except Exception as exc:
                        logging.error('%r generated an exception: %s' % (f_n, exc))
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', ''))
                    file_name += '.sbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True

### multiprocessing

import multiprocessing

## Multiprocessing implementation of the thermodynamics package
#
#
def runThermo_multi_process(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15):
    rpcache = rpCache.rpCache()
    cc_preprocess = rpcache.getCCpreprocess()
    kegg_dG = rpcache.getKEGGdG()
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        with tempfile.TemporaryDirectory() as tmpInputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            ### construct the processes list and start
            processes = []
            for s_l in chunkIt(glob.glob(tmpInputFolder+'/*'), num_workers):
                p = multiprocessing.Process(target=singleThermo, args=(s_l, pathway_id, tmpOutputFolder, kegg_dG, cc_preprocess, ph, ionic_strength, pMg, temp_k))
                processes.append(p)
                p.start()
            #wait for all to finish
            for process in processes:
                process.join()
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    file_name = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', ''))
                    file_name += '.sbml.xml'
                    info = tarfile.TarInfo(file_name)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True
'''


############################# single core ##########################

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

    cc = initThermo()

    rpsbml = rpSBML(inFile)

    results = pathway(rpsbml=rpsbml, cc=cc, pathway_id=pathway_id,
                      ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k,
                      update_rpsbml=True) # ignore the results since written to SBML file
 
    rpsbml.writeSBML(outFile)

    # #mnx_default_conc = json.load(open('data/mnx_default_conc.json', 'r'))
    # mnx_default_conc = json.load(open(os.path.join(os.path.dirname(os.path.abspath( __file__ )), 'data', 'mnx_default_conc.json'), 'r'))


    rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)

    logger.debug('Passing the sbml file: '+str(sbml_path))
    fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '')
    rpequilibrator.rpsbml = rpsbml

    return True


def initThermo():

    cc = ComponentContribution()
    cc.p_h = Q_(ph)
    cc.ionic_strength = Q_(str(ionic_strength)+' mM')
    cc.p_mg = Q_(pMg)
    cc.temperature = Q_(str(temp_k)+' K')

    return cc

def runMDF_hdd(inputTar, outputTar, pathway_id='rp_pathway', thermo_id='dfG_prime_o', fba_id='fba_obj_fraction', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
    """Given a tar input file, perform MDF analysis for each rpSBML file.

    :param inputTar: The path to the input TAR file
    :param outputTar: The path to the output TAR file
    :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
    :param thermo_id: The id of the thermodynamics id (Default: dfG_prime_o)
    :param fba_id: The id of the FBA value (Default: fba_obj_fraction)
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
    :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

    :type inputTar: str
    :type outputTar: str
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
    rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)
    with tempfile.TemporaryDirectory() as tmpInputFolder:
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                logging.debug('=========== '+str(sbml_path)+' ============')
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '') 
                rpsbml = rpSBML.rpSBML(fileName, path=sbml_path)
                rpequilibrator.rpsbml = rpsbml
                res = rpequilibrator.MDF(pathway_id, thermo_id, fba_id, stdev_factor, True) #ignore the results since written to SBML file
                #ignore res since we are passing write to SBML
                rpsbml.writeSBML(tmpOutputFolder)
                rpsbml = None
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', ''))
                    fileName += '.sbml.xml'
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True


def runEqSBtab_hdd(inputTar, outputTar, pathway_id='rp_pathway', fba_id=None, thermo_id='dfG_prime_o', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
    """Given a tar input file, perform MDF analysis for each rpSBML file.

    :param inputTar: The path to the input TAR file
    :param outputTar: The path to the output TAR file
    :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
    :param fba_id: The id of the FBA value. Default sets all FBA values to 1.0 and if specified (Default: None)
    :param thermo_id: The id of the thermodynamics id (Default: dfG_prime_o)
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)
    :param stdev_factor: The standard deviation factor to calculate MDF (Default: 1.96)

    :type inputTar: str
    :type outputTar: str
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
    rpequilibrator = rpEquilibrator.rpEquilibrator(ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k)
    with tempfile.TemporaryDirectory() as tmpInputFolder:
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            tar = tarfile.open(inputTar, mode='r')
            tar.extractall(path=tmpInputFolder)
            tar.close()
            if len(glob.glob(tmpInputFolder+'/*'))==0:
                logging.error('Input file is empty')
                return False
            for sbml_path in glob.glob(tmpInputFolder+'/*'):
                logging.debug('=========== '+str(sbml_path)+' ============')
                fileName = sbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '') 
                rpsbml = rpSBML.rpSBML(fileName, path=sbml_path)
                rpequilibrator.rpsbml = rpsbml
                status = rpequilibrator.toNetworkSBtab(os.path.join(tmpOutputFolder, fileName+'.tsv'), pathway_id, thermo_id, fba_id, stdev_factor)
                rpsbml = None
            if len(glob.glob(tmpOutputFolder+'/*'))==0:
                logging.error('rpThermo has not produced any results')
                return False
            with tarfile.open(outputTar, mode='w:gz') as ot:
                for sbml_path in glob.glob(tmpOutputFolder+'/*'):
                    fileName = str(sbml_path.split('/')[-1])
                    info = tarfile.TarInfo(fileName)
                    info.size = os.path.getsize(sbml_path)
                    ot.addfile(tarinfo=info, fileobj=open(sbml_path, 'rb'))
    return True
    

""" DEPRECATED: only useful when using multiprocessing
##
#
#
#def main(inputTar, outputTar, num_workers=10, pathway_id='rp_pathway', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
def main_thermo(inputTar, outputTar, pathway_id='rp_pathway', ph=7.0, ionic_strength=200, pMg=10.0, temp_k=298.15, stdev_factor=1.96):
    with tempfile.TemporaryDirectory() as tmpCountFolder:
        num_models = 0
        tar = tarfile.open(inputTar, mode='r')
        tar.extractall(path=tmpCountFolder)
        num_models = len(glob.glob(tmpCountFolder+'/*'))
        tar.close()
        if num_workers<=0:
            logging.error('Cannot have less or 0 workers')
            return False
        return runThermo_hdd(inputTar, outputTar, pathway_id=pathway_id, ph=ph, ionic_strength=ionic_strength, pMg=pMg, temp_k=temp_k, stdev_factor=stdev_factor):
        ''' Seems like subprocessing + Equilibrator do not play well together -- needs further testing
        #TODO: count the number of models in the tar and if <num_workers then choose that
        if num_models==0:
            logging.warning('The input tar file seems to be empty') 
            return False
        if num_workers==1 or num_models==1:
            return runThermo_hdd(inputTar, outputTar, pathway_id)
        else:
            #runThermo_multi_concurrent(inputTar, outputTar, num_workers, pathway_id)
            return runThermo_multi_process(inputTar, outputTar, num_workers, pathway_id)
        '''
"""