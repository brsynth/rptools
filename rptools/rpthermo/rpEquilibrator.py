from equilibrator_api                      import Q_
from equilibrator_assets.generate_compound import get_or_create_compound
from equilibrator_assets.group_decompose   import GroupDecompositionError
from equilibrator_pathway                  import Pathway
from equilibrator_cache                    import exceptions as equilibrator_exceptions
import numpy as np
import logging


# WARNING: taking the sum of the reaction thermodynamics is perhaps not the best way to do it
def pathway(rpsbml, cc, pathway_id='rp_pathway', update_rpsbml=True, logger=logging.getLogger(__name__)):
    """Calculate the dG of a heterologous pathway

    :param rpsbml: An rpSBML object
    :param cc: Component contribution
    :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
    :param update_rpsbml: Write the results to the rpSBML file (Default: True)

    :type rpsbml: rptools.rpSBML
    :type cc: ComponentContribution
    :type pathway_id: str
    :type update_rpsbml: bool

    :rtype: tuple
    :return: Tuple with the following information, in order: sum dG_prime, std dG_prime, sum dG_prime_o, std dG_prime_o, sum dG_prime_m, std dG_prime_m. Also False if function error.
    """

    rp_pathway = rpsbml.getModel().getPlugin('groups').getGroup(pathway_id)

    if not rp_pathway:
        logger.error('Cannot retreive the pathway: '+str(pathway_id))
        return False

    pathway_balanced = []
    pathway_reversibility_index = []
    pathway_reversibility_index_error = []
    pathway_standard_dg = []
    pathway_standard_dg_error = []
    pathway_standard_dg_prime = []
    pathway_standard_dg_prime_error = []
    pathway_physiological_dg_prime = []
    pathway_physiological_dg_prime_error = []

    calc_cmp = {}

    for react in [rpsbml.getModel().getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:

        logger.debug('Sending the following reaction to _reactionStrQuery: '+str(react))

        res = _reactionStrQuery(react, rpsbml, cc, update_rpsbml, logger=logger)

        logger.debug('The result is :', res)

        if res != {}:

            # WARNING: the uncertainty for the three thermo calculations should be the same
            pathway_balanced.append(res['rxn_is_balanced'])

            try:
                pathway_reversibility_index.append(res['ln_reversibility_index']['value'])
            except KeyError:
                pathway_reversibility_index.append(0.0)

            try:
                pathway_reversibility_index_error.append(res['ln_reversibility_index']['error'])
            except KeyError:
                pathway_reversibility_index_error.append(0.0)

            # ignoring --  need to see if legacy component contribution can return these values
            # pathway_standard_dg.append(res[2][0])
            # pathway_standard_dg_error.append(res[2][1])
            pathway_standard_dg_prime.append(res['standard_dg_prime']['value'])
            pathway_standard_dg_prime_error.append(res['standard_dg_prime']['error'])
            pathway_physiological_dg_prime.append(res['physiological_dg_prime']['value'])
            pathway_physiological_dg_prime_error.append(res['physiological_dg_prime']['error'])

        else:

            logger.info('Native equilibrator string query failed')
            logger.info('Trying equilibrator_api component contribution')
            logger.debug('Trying to calculate using CC: '+str(react))

            res = _reactionCmpQuery(react, rpsbml, cc, calc_cmp, update_rpsbml, logger=logger)

            if res:
                pathway_standard_dg_prime.append(res[0])
                pathway_standard_dg_prime_error.append(res[2])
                pathway_physiological_dg_prime.append(res[1])
                pathway_physiological_dg_prime_error.append(res[2])
                # TODO: need to implement
                pathway_balanced.append(None)
                pathway_reversibility_index.append(None)
                pathway_reversibility_index_error.append(None)
            else:
                logger.warning('Cannot calculate the thermodynmics for the reaction: '+str(react))
                logger.warning('Setting everything to 0')
                pathway_standard_dg_prime.append(0.0)
                pathway_standard_dg_prime_error.append(0.0)
                pathway_physiological_dg_prime.append(0.0)
                pathway_physiological_dg_prime_error.append(0.0)
                # TODO: need to implement
                pathway_balanced.append(None)
                pathway_reversibility_index.append(None)
                pathway_reversibility_index_error.append(None)
                return {}

    # WARNING return is ignoring balanced and reversibility index -- need to implement in legacy to return it (however still writing these results to the SBML)
    if update_rpsbml:
        rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_o'     , np.sum(pathway_standard_dg_prime)        , 'kj_per_mol')
        rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_o_std' , np.std(pathway_standard_dg_prime)        , 'kj_per_mol')
        rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_m'     , np.sum(pathway_physiological_dg_prime)   , 'kj_per_mol')
        rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_prime_m_std' , np.std(pathway_physiological_dg_prime)   , 'kj_per_mol')
        rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_uncert'      , np.mean(pathway_standard_dg_prime_error) , 'kj_per_mol')
        rpsbml.addUpdateBRSynth(rp_pathway, 'dfG_uncert_std'  , np.std(pathway_standard_dg_prime_error)  , 'kj_per_mol')

    return {
        'dfG_prime_o'     : np.sum(pathway_standard_dg_prime),
        'dfG_prime_o_std' : np.std(pathway_standard_dg_prime),
        'dfG_prime_m'     : np.sum(pathway_physiological_dg_prime),
        'dfG_prime_m_std' : np.std(pathway_physiological_dg_prime),
        'dfG_uncert'      : np.mean(pathway_standard_dg_prime_error),
        'dfG_uncert_std'  : np.std(pathway_standard_dg_prime_error)
    }


# TODO: when an inchikey is passed, (and you don't have any other xref) and equilibrator finds the correct species then update the MIRIAM annotations
def _reactionStrQuery(libsbml_reaction, rpsbml, cc, update_rpsbml=False, logger=logging.getLogger(__name__)):
    """Build the string reaction from a libSBML reaction object to send to equilibrator and return the different thermodynamics analysis available

    :param libsbml_reaction: A libsbml reaction object
    :param rpsbml: An rpSBML object
    :param cc: Component contribution
    :param write_results: Write the results to the rpSBML file (Default: False)

    :type libsbml_reaction: libsbml.Reaction
    :type rpsbml: rptools.rpSBML
    :type cc: ComponentContribution
    :type write_results: bool

    :rtype: bool
    :return: Success or failue of the function
    """

    reac_str = ''

    results = {}

    try:
        reac_str = _makeReactionStr(libsbml_reaction, rpsbml, logger=logger)
        logger.debug('The reaction string is: '+str(reac_str))
        if not reac_str:
            logger.warning('Could not generate the reaction string for: '+str(libsbml_reaction))
            if update_rpsbml:
                logger.warning('Writing the 0 results to the file')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_o', 0.0, 'kj_per_mol')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_m', 0.0, 'kj_per_mol')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_uncert', 0.0, 'kj_per_mol')
                #Are there default values for these?
                #rpsbml.addUpdateBRSynth(libsbml_reaction, 'reversibility_index', 0.0)
                #rpsbml.addUpdateBRSynth(libsbml_reaction, 'balanced', rxn.is_balanced())
            return {}

        rxn = cc.parse_reaction_formula(reac_str)
        results['standard_dg']            = cc.standard_dg(rxn)
        results['standard_dg_prime']      = cc.standard_dg_prime(rxn)
        results['physiological_dg_prime'] = cc.physiological_dg_prime(rxn)
        results['ln_reversibility_index'] = cc.ln_reversibility_index(rxn)

        if type(ln_reversibility_index) == float:
            logger.warning('The reversibility index is infinite: '+str(ln_reversibility_index))
            results['ln_reversibility_index']       = None
            results['ln_reversibility_index_error'] = None

        else:
            ln_reversibility_index_error = ln_reversibility_index.error.m
            ln_reversibility_index = ln_reversibility_index.value.m

        logger.debug(rxn.is_balanced())
        logger.debug('ln_reversibility_index:         ' + str(results['ln_reversibility_index']))
        logger.debug('standard_dg.value.m:            ' + str(results['standard_dg'].value.m))
        logger.debug('standard_dg.error.m:            ' + str(results['standard_dg'].error.m))
        logger.debug('standard_dg_prime.value.m:      ' + str(results['standard_dg_prime'].value.m))
        logger.debug('standard_dg_prime.error.m:      ' + str(results['standard_dg_prime'].error.m))
        logger.debug('physiological_dg_prime.value.m: ' + str(results['physiological_dg_prime'].value.m))
        logger.debug('physiological_dg_prime.error.m: ' + str(results['physiological_dg_prime'].error.m))

        if update_rpsbml:
            rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_o' , results['standard_dg'].value.m            , 'kj_per_mol')
            rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_m' , results['physiological_dg_prime'].value.m , 'kj_per_mol')
            rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_uncert'  , standard_dg.error.m                       , 'kj_per_mol')
            if results['ln_reversibility_index']:
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'reversibility_index', results['ln_reversibility_index'])
            rpsbml.addUpdateBRSynth(libsbml_reaction, 'balanced', rxn.is_balanced())

        return {
            'rxn_is_balanced': rxn.is_balanced(),
            'ln_reversibility_index': {
                'value': results['ln_reversibility_index'],
                'error': results['ln_reversibility_index_error']
                },
            'standard_dg': {
                'value': float(results['standard_dg'].value.m),
                'error': float(results['standard_dg'].error.m)
                },
            'standard_dg_prime': {
                'value': float(results['standard_dg_prime'].value.m),
                'error': float(results['standard_dg_prime'].error.m)
                },
            'physiological_dg_prime': {
                'value': float(results['physiological_dg_prime'].value.m),
                'error': float(results['physiological_dg_prime'].error.m)
            }
        }

        # return (rxn.is_balanced(),
        #         (ln_reversibility_index, ln_reversibility_index_error),
        #         (float(standard_dg.value.m), float(standard_dg.error.m)),
        #         (float(standard_dg_prime.value.m), float(standard_dg_prime.error.m)), 
        #         (float(physiological_dg_prime.value.m), float(physiological_dg_prime.error.m)))
    except equilibrator_exceptions.ParseException:
        logger.warning('One of the reaction species cannot be parsed by equilibrator: '+str(reac_str))
        return {}

    except equilibrator_exceptions.MissingDissociationConstantsException:
        logger.warning('Some of the species have not been pre-caclulated using ChemAxon')



# # TODO: need to report when calculating the thermodynamics of reactions failed.... perhaps in the pathway add True/False tag to see
# class rpEquilibrator:
#     """Class containing collection of functions to intereact between rpSBML files and equilibrator. Includes a function to convert an rpSBML file to a SBtab format for MDF analysis
#     """
#     def __init__(self, rpsbml=None, ph=7.5, ionic_strength=200, pMg=10.0, temp_k=298.15):
#         """Constructor class for rpEquilibrator

#         :param rpsbml: rpSBML object (Default: None)
#         :param ph: pH of the cell input from the rpSBML model input (Default: 7.5)
#         :param ionic_strength: Ionic strength from the rpSBML model input (Default: 200)
#         :param pMg: pMg value from the rpSBML model input (Default: 10.0)
#         :param temp_k: Temperature from the rpSBML model input in Kelvin (Default: 298.15)

#         :type rpsbml: rpSBML
#         :type ph: float
#         :type ionic_strength: float
#         :type pMg: float
#         :type temp_k: float

#         .. document private functions
#         .. automethod:: _makeSpeciesStr
#         .. automethod:: _makeReactionStr
#         .. automethod:: _speciesCmpQuery
#         .. automethod:: _reactionCmpQuery
#         .. automethod:: _reactionStrQuery
#         """
#         logger = logging.getLogger(__name__)
#         logger.debug('Started instance of rpEquilibrator')
#         cc = ComponentContribution()
#         cc.p_h = Q_(ph)
#         cc.ionic_strength = Q_(str(ionic_strength)+' mM')
#         cc.p_mg = Q_(pMg)
#         cc.temperature = Q_(str(temp_k)+' K')
#         ph = ph
#         ionic_strength = ionic_strength
#         pMg = pMg
#         temp_k = temp_k
#         #mnx_default_conc = json.load(open('data/mnx_default_conc.json', 'r'))
#         mnx_default_conc = json.load(open(os.path.join(os.path.dirname(os.path.abspath( __file__ )), 'data', 'mnx_default_conc.json'), 'r'))
#         rpsbml = rpsbml
#         calc_cmp = {}
    

#     ##################################################################################
#     ############################### PRIVATE ##########################################
#     ##################################################################################


# TODO: metanetx.chemical:MNXM7 + bigg.metabolite:pi
def _makeSpeciesStr(libsbml_species, rpsbml, ret_type='xref', logger=logging.getLogger(__name__)):
    """Private function that makes a Equilibrator friendly string of a species

    :param libsbml_species: A libsbml species object
    :param ret_type: Type of output. Valid output include: ['name', 'id', 'xref']

    :type libsbml_species: libsbml.Species
    :type ret_type: str

    Take a libsbml species object, parse the MIRIAM or the brsynth (if present) to return 
    the equilibrator appropriate string. The order of preference is the following:
    example input MIRIAM annotation: {'inchikey': ['GPRLSGONYQIRFK-UHFFFAOYSA-N'], 'seed': ['cpd00067'], 'sabiork': ['39'], 'reactome': ['R-ALL-74722', 'R-ALL-70106', 'R-ALL-5668577', 'R-ALL-428548', 'R-ALL-428040', 'R-ALL-427899', 'R-ALL-425999', 'R-ALL-425978', 'R-ALL-425969', 'R-ALL-374900', 'R-ALL-372511', 'R-ALL-351626', 'R-ALL-2872447', 'R-ALL-2000349', 'R-ALL-194688', 'R-ALL-193465', 'R-ALL-163953', 'R-ALL-156540', 'R-ALL-1470067', 'R-ALL-113529', 'R-ALL-1132304'], 'metacyc': ['PROTON'], 'hmdb': ['HMDB59597'], 'chebi': ['5584', '13357', '10744', '15378'], 'bigg': ['M_h', 'h'], 'metanetx': ['MNXM89553', 'MNXM145872', 'MNXM1', 'MNXM01']}
    -KEGG
    -CHEBI
    #-bigg
    #-MNX
    -inchikey
    ret_type -> valid options (xref, id, name)

    :rtype: str
    :return: The string id of the species or False if fail
    """

    logger.debug('ret_type: '+str(ret_type))

    if ret_type=='name':
        return libsbml_species.getName()

    elif ret_type=='id':
        return libsbml_species.getId()

    elif ret_type=='xref':
        annot = libsbml_species.getAnnotation()
        if not annot:
            logger.error('Cannot retreive the annotation')
            return False
        miriam_dict = rpsbml.readMIRIAMAnnotation(annot)
        logger.debug('miriam_dict: '+str(miriam_dict))
        if not miriam_dict:
            logger.warning('The object annotation does not have any MIRIAM entries')
            return False
        if 'kegg' in miriam_dict:
            if miriam_dict['kegg']:
                try:
                    #take the lowest value
                    int_list = [int(i.replace('C', '')) for i in miriam_dict['kegg']]
                    return 'KEGG:'+str(miriam_dict['kegg'][int_list.index(min(int_list))])
                except ValueError:
                    logger.warning('There is a non int value in: '+str(miriam_dict['kegg']))
        elif 'chebi' in miriam_dict:
            if miriam_dict['chebi']:
                try:
                    #take the lowest value
                    int_list = [int(i) for i in miriam_dict['chebi']]
                    return 'CHEBI:'+str(miriam_dict['chebi'][int_list.index(min(int_list))])
                except ValueError:
                    logger.warning('There is a non int value in: '+str(miriam_dict['chebi']))
        elif 'metanetx' in miriam_dict:
            if miriam_dict['metanetx']:
                try:
                    #take the lowest value
                    int_list = [int(i.replace('MNXM', '')) for i in miriam_dict['metanetx']]
                    return 'metanetx.chemical:'+str(miriam_dict['metanetx'][int_list.index(min(int_list))])
                except ValueError:
                    logger.warning('There is a non int value in: '+str(miriam_dict['metanetx']))
        elif 'inchikey' in miriam_dict:
            if miriam_dict['inchikey']:
                if len(miriam_dict['inchikey'])==1:
                    return miriam_dict['inchikey'][0]
                else:
                    logger.warning('There are multiple values of inchikey: '+str(miriam_dict['inchikey']))
                    logger.warning('Taking the first one')
                    return miriam_dict['inchikey'][0]
        else:
            logger.warning('Could not extract string input for '+str(miriam_dict))
            return False
        logger.warning('The MIRIAM annotation does not have the required information')
        return False
    else:
        logger.warning('Cannot determine ret_type: '+str(ret_type))


def _makeReactionStr(libsbml_reaction, rpsbml, ret_type='xref', ret_stoichio=True, logger=logging.getLogger(__name__)):
    """Make the reaction formulae string to query equilibrator

    :param libsbml_reaction: A libsbml reaction object
    :param rpsbml: An rpSBML object
    :param ret_type: Type of output. Valid output include: ['name', 'id', 'xref'] (Default: xref)
    :param ret_stoichio: Return the stoichio or not (Default: True)

    :type libsbml_reaction: libsbml.Reaction
    :type rpsbml: rptools.rpSBML
    :type ret_type: str
    :type ret_stoichio: bool

    :rtype: str
    :return: The string id of the reaction or False if fail
    """

    reac_str = ''

    for rea in libsbml_reaction.getListOfReactants():
        rea_str = _makeSpeciesStr(rpsbml.getModel().getSpecies(rea.getSpecies()), rpsbml, ret_type, logger=logger)
        if rea_str:
            if ret_stoichio:
                reac_str += str(rea.getStoichiometry())+' '+str(rea_str)+' + '
            else:
                reac_str += str(rea_str)+' + '
        else:
            return False

    reac_str = reac_str[:-2]
    reac_str += '<=> ' # TODO: need to find a way to determine the reversibility of the reaction

    for pro in libsbml_reaction.getListOfProducts():
        pro_str = _makeSpeciesStr(rpsbml.getModel().getSpecies(pro.getSpecies()), rpsbml, ret_type, logger=logger)
        if pro_str:
            if ret_stoichio:
                reac_str += str(pro.getStoichiometry())+' '+str(pro_str)+' + '
            else:
                reac_str += str(pro_str)+' + '
        else:
            return False

    reac_str = reac_str[:-2]

    logger.debug('reac_str: '+str(reac_str))

    return reac_str


################### Equilibrator component contribution queries instead of using the native functions ###########


def _speciesCmpQuery(libsbml_species, rpsbml, cc, calc_cmp, logger=logging.getLogger(__name__)):
    """Use the native equilibrator-api compound contribution method

    :param libsbml_species: A libsbml species object
    :param rpsbml: An rpSBML object
    :param cc: Component contribution

    :type libsbml_species: libsbml.Reaction
    :type rpsbml: rptools.rpSBML
    :type cc: ComponentContribution

    :rtype: tuple
    :return: Tuple of size two with mu and sigma values in that order or (None, None) if fail
    """

    annot = libsbml_species.getAnnotation()

    if not annot:
        logger.warning('The annotation of '+str(libsbml_species)+' is None....')
        return None, None

    brs_annot = rpsbml.readBRSYNTHAnnotation(libsbml_species.getAnnotation())

    # TODO: handle the condition where there are no inchi values but there are SMILES -- should rarely, if ever happen
    logger.debug('libsbml_species: '+str(libsbml_species))
    # logger.debug('brs_annot: '+str(brs_annot))
    # Try to get the cmp from the ID
    spe_id = _makeSpeciesStr(libsbml_species, rpsbml, logger=logger)
    spe_cmp = None

    if spe_id:
        logger.debug('Trying to find the CMP using the xref string: '+str(spe_id))
        spe_cmp = cc.ccache.get_compound(_makeSpeciesStr(libsbml_species, rpsbml, logger=logger))

    if not spe_cmp:
        logger.debug('Trying to find the CMP using the structure')
        # try to find it in the local data - we do this because there are repeated species in many files
        if brs_annot['inchi'] in calc_cmp:
            spe_cmp = calc_cmp[brs_annot['inchi']]
        elif brs_annot['smiles'] in calc_cmp:
            spe_cmp = calc_cmp[brs_annot['smiles']]
        else:
            # if you cannot find it then calculate it
            try:
                spe_cmp = get_or_create_compound(cc.ccache, brs_annot['inchi'], mol_format='inchi')
                calc_cmp[brs_annot['inchi']] = spe_cmp
                calc_cmp[brs_annot['smiles']] = spe_cmp
            except (OSError, KeyError, GroupDecompositionError) as e:
                try:
                    spe_cmp = get_or_create_compound(cc.ccache, brs_annot['smiles'], mol_format='smiles')
                    calc_cmp[brs_annot['smiles']] = spe_cmp
                    calc_cmp[brs_annot['inchi']] = spe_cmp
                except (OSError, KeyError, GroupDecompositionError) as e:
                    logger.warning('The following species does not have brsynth annotation InChI or SMILES: '+str(libsbml_species.getId()))
                    logger.warning('Or Equilibrator could not convert the structures')
                    logger.warning(e)
                    return None, None

    if spe_cmp.id == 4: # this is H+ and can be ignored
        return 'h', 'h'

    logger.debug('spe_cmp: '+str(spe_cmp))

    # mu -- The first value is the formation energy mean estimate (as before)
    # uncertainty_array -- The second value is an array representing the uncertainty (a projection in the subspace of known reactions). If you take the inner product of that vector with itself, you'll get the variance of the estimate. An inner product with the array of another compound would yield the co-variance of their estimates.
    # null_space -- The third return-value is similar to the second, but is a projection on the null-space. If the inner product with another such vector is not zero, it means the estimate is completely unreliable and should not be used.
    mu, uncertainty_array, null_space = cc.predictor.preprocess.get_compound_prediction(spe_cmp)

    return mu, math.sqrt(numpy.inner(uncertainty_array))


def _reactionCmpQuery(libsbml_reaction, rpsbml, cc, calc_cmp, update_rpsbml=False, physio_param=1e-3, logger=logging.getLogger(__name__)):
    """This method makes a list of structure compounds and uses equilibrator to return the reaction dG

    :param libsbml_reaction: A libsbml reaction object
    :type rpsbml: rptools.rpSBML
    :param cc: Component contribution
    :param write_results: Write the results to the rpSBML file (Default: False)
    :param physio_param: The physiological parameter, i.e. the concentration of the compounds to calculate the dG (Default: 1e-3)

    :type libsbml_reaction: libsbml.Reaction
    :param rpsbml: An rpSBML object
    :type cc: ComponentContribution
    :type write_results: bool
    :type physio_param: float

    :rtype: tuple
    :return: Tuple of size three with dfG_prime_o, dfG_prime_m, uncertainty values in that order or False if fail
    """

    mus = []
    sigma_vecs = []
    S = []
    dfG_prime_o = None
    dfG_prime_m = None
    uncertainty = None

    for rea in libsbml_reaction.getListOfReactants():
        logger.debug('------------------- '+str(rea.getSpecies())+' --------------')
        mu, sigma = _speciesCmpQuery(rpsbml.getModel().getSpecies(rea.getSpecies()), rpsbml, cc, calc_cmp, logger=logger)
        logger.debug('mu: '+str(mu))
        if not mu:
            logger.warning('Failed to calculate the reaction mu thermodynamics using compound query')
            if update_rpsbml:
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_o', 0.0, 'kj_per_mol')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_m', 0.0, 'kj_per_mol')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_uncert', 0.0, 'kj_per_mol')
            return False
        elif mu == 'h': # skipping the Hydrogen
            continue
        mus.append(mu)
        sigma_vecs.append(sigma)
        S.append([-rea.getStoichiometry()])

    for pro in libsbml_reaction.getListOfProducts():
        mu, sigma = _speciesCmpQuery(rpsbml.getModel().getSpecies(pro.getSpecies()), rpsbml, cc, calc_cmp, logger=logger)
        if not mu:
            logger.warning('Failed to calculate the reaction mu thermodynamics using compound query')
            if update_rpsbml:
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_o', 0.0, 'kj_per_mol')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_m', 0.0, 'kj_per_mol')
                rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_uncert', 0.0, 'kj_per_mol')
            return False
        elif mu == 'h': # skipping the Hydrogen
            continue
        mus.append(mu)
        sigma_vecs.append(sigma)
        S.append([pro.getStoichiometry()])

    mus = Q_(mus, 'kJ/mol')
    sigma_vecs = Q_(sigma_vecs, 'kJ/mol')
    np_S = np.array(S)
    dfG_prime_o = np_S.T@mus
    dfG_prime_o = float(dfG_prime_o.m[0])
    ###### adjust fot physio parameters to calculate the dGm'
    # TODO: check with Elad
    dfG_prime_m = float(dfG_prime_o)+float(cc.RT.m)*sum([float(sto[0])*float(np.log(co)) for sto, co in zip(S, [physio_param]*len(S))])
    uncertainty = np_S.T@sigma_vecs
    uncertainty = uncertainty@uncertainty.T
    uncertainty = uncertainty.m[0][0]

    if update_rpsbml:
        rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_o', dfG_prime_o, 'kj_per_mol')
        rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_prime_m', dfG_prime_m, 'kj_per_mol')
        rpsbml.addUpdateBRSynth(libsbml_reaction, 'dfG_uncert', uncertainty, 'kj_per_mol')

    return dfG_prime_o, dfG_prime_m, uncertainty


'''
## Not sure if we should implement such a function -- recommended by Elad I geuss
#
#
def pathwayCmpQuery(self, write_results=False):
    #1) build the stochio matrix taking into account all the species of the reaction -- must keep track
    pass

#################### native equilibrator-api functions ###############

def speciesStrQuery(self, libsbml_species, write_results=False):
    """
    Return the formation energy of a chemical species
    """
    return False
'''

################################################################################
########################### PUBLIC FUNCTIONS ###################################
################################################################################




def toNetworkSBtab(self, output, pathway_id='rp_pathway', thermo_id='dfG_prime_o', fba_id='fba_obj_fraction', stdev_factor=1.96):
    """Convert an SBML pathway to a simple network for input to equilibrator-pathway for MDF

    :param output: Output path of the TSV file
    :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
    :param thermo_id: The id of the thermodynamics result to be exported to the SBtab file (Default: dfG_prime_o, Valid Options: [dfG_prime_o, dfG_prime_m])
    :param fba_id: The id of the FBA value to be exported to SBtab (Default: fba_obj_fraction)
    :param stdev_factor: The standard deviation factor (Default: 1.96)

    :type output: str
    :type pathway_id: str
    :type thermo_id: str
    :type fba_id: str
    :type stdev_factor: float

    :rtype: bool
    :return: Success or failure of the function
    """
    groups = rpsbml.getModel().getPlugin('groups')
    rp_pathway = groups.getGroup(pathway_id)
    if not rp_pathway:
        logger.error('Cannot retreive the pathway: '+str(pathway_id))
        return False
    with open(output, 'w') as fo:
        ####################### Make the header of the document ##############
        fo.write("!!!SBtab DocumentName='E. coli central carbon metabolism - balanced parameters' SBtabVersion='1.0'\t\t\t\n")
        fo.write("!!SBtab TableID='Configuration' TableType='Config'\t\t\t\n")
        fo.write("!Option\t!Value\t!Comment\t\n")
        fo.write("algorithm\tMDF\tECM, or MDF\t\n")
        fo.write("p_h\t"+str(ph)+"\t\t\n")
        fo.write("ionic_strength\t"+str(ionic_strength)+" mM\t\t\n")
        fo.write("p_mg\t"+str(pMg)+"\t\t\n")
        fo.write("stdev_factor    "+str(stdev_factor)+"\n")
        fo.write("\t\t\t\n")
        ####################### Make the reaction list ######################
        fo.write("!!SBtab TableID='Reaction' TableType='Reaction'\t\t\t\n")
        fo.write("!ID\t!ReactionFormula\t\t\n")
        #TODO: need to sort it in the reverse order
        #TODO: use the rpGraph to sort the pathway in the right order
        ordered_react = [rpsbml.getModel().getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]
        ordered_react.sort(key=lambda x: int(x.getId().replace('RP', '')), reverse=True)
        #for react in [rpsbml.getModel().getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:   
        for react in ordered_react:
            react_str = _makeReactionStr(react, 'id', True)
            if react_str:
                fo.write(str(react.getId())+"\t"+str(react_str)+"\n")
            else:
                logger.error('Cannot build the reaction: '+str(rect))
                return False
        fo.write("\t\t\t\n")
        fo.write("\t\t\t\n")
        ########################## Make the species list ###################
        fo.write("!!SBtab TableID='Compound' TableType='Compound'\t\t\t\n")
        fo.write("!ID\t!Identifiers\t\t\n")
        rp_species = rpsbml.readUniqueRPspecies(pathway_id)
        for spe_id in rp_species:
            spe = rpsbml.getModel().getSpecies(spe_id)
            miriam_dict = rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
            if not miriam_dict:
                logger.warning('The object annotation does not have any MIRIAM entries')
                return False
            iden_str = None
            if 'kegg' in miriam_dict:
                if miriam_dict['kegg']:
                    try:
                        #take the lowest value
                        int_list = [int(i.replace('C', '')) for i in miriam_dict['kegg']]
                        iden_str = 'KEGG:'+str(miriam_dict['kegg'][int_list.index(min(int_list))])
                    except ValueError:
                        logger.warning('There is a non int value in: '+str(miriam_dict['kegg']))
            if 'chebi' in miriam_dict and not iden_str:
                if miriam_dict['chebi']:
                    try:
                        #take the lowest value
                        int_list = [int(i) for i in miriam_dict['chebi']]
                        iden_str = 'CHEBI:'+str(miriam_dict['chebi'][int_list.index(min(int_list))])
                    except ValueError:
                        logger.warning('There is a non int value in: '+str(miriam_dict['chebi']))
            if 'metanetx' in miriam_dict and not iden_str:
                if miriam_dict['metanetx']:
                    try:
                        #take the lowest value
                        int_list = [int(i.replace('MNXM', '')) for i in miriam_dict['metanetx']]
                        iden_str = 'metanetx.chemical:'+str(miriam_dict['metanetx'][int_list.index(min(int_list))])
                    except ValueError:
                        logger.warning('There is a non int value in: '+str(miriam_dict['metanetx']))
            if 'inchikey' in miriam_dict and not iden_str:
                if miriam_dict['inchikey']:
                    if len(miriam_dict['inchikey'])==1:
                        iden_str = miriam_dict['inchikey'][0]
                    else:
                        logger.warning('There are multiple values of inchikey: '+str(miriam_dict['inchikey']))
                        logger.warning('Taking the first one')
                        iden_str = miriam_dict['inchikey'][0]
            if not iden_str:
                logger.warning('Could not extract string input for '+str(miriam_dict))
            fo.write(str(spe_id)+"\t"+str(iden_str)+"\t\t\n")
        fo.write("\t\t\t\n")
        ################## Add FBA values ##############################
        #TODO: perhaps find a better way than just setting this to 1
        fo.write("!!SBtab TableID='Flux' TableType='Quantity' Unit='mM/s'\t\t\t\n")
        fo.write("!QuantityType\t!Reaction\t!Value\t\n")
        for react in [rpsbml.getModel().getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:
            brs_annot = rpsbml.readBRSYNTHAnnotation(react.getAnnotation())
            if fba_id:
                if fba_id in brs_annot:
                    #the saved value is mmol_per_gDW_per_hr while rpEq required mmol_per_s
                    #WARNING: rpEq says that these values are mMol/s while here we have mMol/gDw/h. However not changing since this would mean doing /3600 and
                    #given the values that seems wrong
                    fo.write("rate of reaction\t"+str(react.getId())+"\t"+str(brs_annot[fba_id]['value'])+"\t\n")
                else:
                    logger.warning('Cannot retreive the FBA value '+str(fba_id)+'. Setting a default value of 1.')
                    fo.write("rate of reaction\t"+str(react.getId())+"\t1\t\n")
            else:
                fo.write("rate of reaction\t"+str(react.getId())+"\t1\t\n")
        ################## Add the concentration bounds ##############################
        fo.write("\t\t\t\n")
        fo.write("!!SBtab TableID='ConcentrationConstraint' TableType='Quantity' Unit='mM'\t\t\t\n")
        fo.write("!QuantityType\t!Compound\t!Min\t!Max\n")
        for spe_id in rp_species: 
            logger.debug('========= '+str(spe_id)+' ========')
            is_found = False
            spe = rpsbml.getModel().getSpecies(spe_id)
            miriam_dict = rpsbml.readMIRIAMAnnotation(spe.getAnnotation())
            logger.debug(miriam_dict)
            if not miriam_dict:
                logger.warning('The object annotation does not have any MIRIAM entries')
                continue
            if 'metanetx' in miriam_dict:
                logger.debug(miriam_dict['metanetx'])
                for mnx in miriam_dict['metanetx']:
                    if mnx in list(mnx_default_conc.keys()) and not is_found:
                        logger.debug('Found default concentration range for '+str(spe.getId())+' ('+str(mnx)+'): '+str(mnx_default_conc[mnx]))
                        fo.write("concentration\t"+spe.getId()+"\t"+str(mnx_default_conc[mnx]['c_min'])+"\t"+str(mnx_default_conc[mnx]['c_max'])+"\n")
                        is_found = True
            if not is_found:
                logger.debug('Using default values for '+str(spe.getId()))
                fo.write("concentration\t"+spe.getId()+"\t0.001\t10\n")
        fo.write("\t\t\t\n")
        ############################ Add the thermo value ###########################
        #TODO: perform on the fly thermodynamic calculations when the values are not included within the SBML file
        if thermo_id:
            fo.write("!!SBtab TableID='Thermodynamics' TableType='Quantity' StandardConcentration='M'\t\t\t\n")
            fo.write("!QuantityType\t!Reaction\t!Compound\t!Value\t!Unit\n")
            for react in [rpsbml.getModel().getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]:
                brs_annot = rpsbml.readBRSYNTHAnnotation(react.getAnnotation())
                try:
                    #TODO: switch to dfG_prime_m when you are sure how to calculate it using the native equilibrator function
                    if thermo_id in brs_annot:
                        if brs_annot[thermo_id]:
                            fo.write("reaction gibbs energy\t"+str(react.getId())+"\t\t"+str(brs_annot['dfG_prime_o']['value'])+"\tkJ/mol\n")
                        else:
                            logger.error(str(thermo_id)+' is empty. Was rpThermodynamics run on this SBML? Aborting...')
                            return False
                    else:
                        logger.error('There is no '+str(thermo_id)+' in the reaction '+str(react.getId()))
                        return False
                except KeyError:
                    logger.error('The reaction '+str(react.getId())+' does not seem to have the following thermodynamic value: '+str(thermo_id))
                    return False
    return True


def MDF(self, pathway_id='rp_pathway', thermo_id='dfG_prime_o', fba_id='fba_obj_fraction', stdev_factor=1.96, write_results=True):
    """Perform MDF analysis on the rpSBML file 

    :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
    :param thermo_id: The id of the thermodynamics result to be exported to the SBtab file (Default: dfG_prime_o, Valid Options: [dfG_prime_o, dfG_prime_m])
    :param fba_id: The id of the FBA value to be exported to SBtab (Default: fba_obj_fraction)
    :param stdev_factor: The standard deviation factor (Default: 1.96)
    :param write_results: Write the results to the rpSBML file (Default: True)

    :type pathway_id: str
    :type thermo_id: str
    :type fba_id: str
    :type stdev_factor: float
    :type write_results: bool

    :rtype: float
    :return: MDF of the pathway
    """
    to_ret_mdf = None
    groups = rpsbml.getModel().getPlugin('groups')
    rp_pathway = groups.getGroup(pathway_id)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        path_sbtab = os.path.join(tmpOutputFolder, 'tmp_sbtab.tsv')
        sbtab_status = toNetworkSBtab(path_sbtab, pathway_id=pathway_id, thermo_id=thermo_id, fba_id=fba_id, stdev_factor=stdev_factor)
        if not sbtab_status:
            logger.error('There was a problem generating the SBtab... aborting')
            return 0.0
        try:
            pp = Pathway.from_sbtab(path_sbtab, comp_contrib=cc)
            pp.update_standard_dgs()
            try:
                mdf_sol = pp.calc_mdf()
            except:
                logger.warning('The calc_mdf function failed')
                logger.warning('Exception: Cannot solve MDF primal optimization problem')
                rpsbml.addUpdateBRSynth(rp_pathway, 'MDF', 0.0, 'kj_per_mol')
                return 0.0
            #mdf_sol = pp.mdf_analysis()
            #plt_reac_plot = mdf_sol.reaction_plot
            #plt_cmp_plot = mdf_sol.compound_plot
            to_ret_mdf = float(mdf_sol.mdf.m)
            if write_results:
                rpsbml.addUpdateBRSynth(rp_pathway, 'MDF', float(mdf_sol.mdf.m), 'kj_per_mol')
        except KeyError as e:
            logger.warning('Cannot calculate MDF')
            logger.warning(e)
            rpsbml.addUpdateBRSynth(rp_pathway, 'MDF', 0.0, 'kj_per_mol')
            return 0.0
        except equilibrator_exceptions.MissingDissociationConstantsException as e:
            logger.warning('Some species are invalid: '+str(e))
            rpsbml.addUpdateBRSynth(rp_pathway, 'MDF', 0.0, 'kj_per_mol')
            return 0.0
    return to_ret_mdf