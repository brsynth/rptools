#!/usr/bin/env python3

from rr_cache import rrCache
from rptools.rplibs import rpSBML
from libsbml        import writeSBMLToFile
import logging


class inchikeyMIRIAM:
    """This class holds the cache information used by the scripts
    """
    def __init__(self, logger=None):
        """Constructor for the inchikeyMIRIAM class
        """
        if logger is None:
            # Create logger
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(getattr(logging, 'ERROR'))
            self.logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')
        else:
            self.logger = logger

        self.logger.info('Started instance of inchikeyMIRIAM')

        self.cache           = rrCache(['deprecatedCID_cid', 'cid_strc', 'chebi_cid'], logger=self.logger)
        self.deprecatedCID_cid = self.cache.get('deprecatedCID_cid')
        self.cid_strc          = self.cache.get('cid_strc')
        self.chebi_cid         = self.cache.get('chebi_cid')


    # def _checkCIDdeprecated(self, cid):
    #     """Function to create return the uniform compound ID
    #
    #     :param cid: The compound id
    #
    #     :type cid: str
    #
    #     :rtype: str
    #     :return: A valid compound id
    #     """
    #     try:
    #         return self.deprecatedCID_cid[cid]
    #     except KeyError:
    #         return cid


    def addInChiKey(self, input_sbml, output_sbml):
        """Check the MIRIAM annotation for MetaNetX or CHEBI id's and try to recover the inchikey from cache and add it to MIRIAM

        :param input_sbml: SBML file input
        :param output_sbml: Output SBML file

        :type input_sbml: str
        :type output_sbml: str

        :rtype: bool
        :return: Success or failure of the function
        """
        filename = input_sbml.split('/')[-1].replace('.rpsbml', '').replace('.sbml', '').replace('.xml', '')
        self.logger.debug(filename)
        rpsbml = rpSBML(inFile=input_sbml, logger=self.logger)
        for spe in rpsbml.getModel().getListOfSpecies():
            inchikey = None
            miriam_dict = rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
            if 'inchikey' in miriam_dict:
                self.logger.info('The species '+str(spe.id)+' already has an inchikey... skipping')
                continue
            try:
                for mnx in miriam_dict['metanetx']:
                    inchikey = self.cid_strc[self.cache._checkCIDdeprecated(mnx, self.deprecatedCID_cid)]['inchikey']
                    if inchikey:
                        rpsbml.addUpdateMIRIAM(spe, 'species', {'inchikey': [inchikey]})
                    else:
                        self.logger.warning('The inchikey is empty for: '+str(spe.id))
                    continue
            except KeyError:
                try:
                    for chebi in miriam_dict['chebi']:
                        inchikey = self.cid_strc[self.cache._checkCIDdeprecated(self.chebi_cid[chebi], self.deprecatedCID_cid)]['inchikey']
                        if inchikey:
                            rpsbml.addUpdateMIRIAM(spe, 'species', {'inchikey': [inchikey]})
                        else:
                            self.logger.warning('The inchikey is empty for: '+str(spe.id))
                        continue
                except KeyError:
                    self.logger.warning('Cannot find the inchikey for: '+str(spe.id))
        writeSBMLToFile(rpsbml.document, output_sbml)
        return True

# def main(input_sbml, output_sbml):
#     """Main function that creates a inchikeyMIRIAM object and runs it
#
#     :param input_sbml: SBML file input
#     :param output_sbml: Output SBML file
#
#     :type input_sbml: str
#     :type output_sbml: str
#
#     :rtype: None
#     :return: None
#     """
#     inchikeymiriam = inchikeyMIRIAM()
#     inchikeymiriam.addInChiKey(input_sbml, output_sbml)
