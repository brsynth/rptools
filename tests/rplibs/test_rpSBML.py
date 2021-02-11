"""
Created on June 17 2020

@author: Joan Hérisson
"""

from    rptools.rplibs import rpSBML
from          tempfile import (
    NamedTemporaryFile
)
from           pathlib import Path
from                os import path  as os_path
from              json import load  as json_load
from              json import loads as json_loads
from                io import open  as io_open
from main import Main_rplibs
from brs_utils         import extract_gz


class Test_rpSBML(Main_rplibs):


    # To avoid limit in dictionaries comparison
    maxDiff = None
    gem_path_gz = os_path.join(
        Main_rplibs.data_path,
        'gem.xml.gz'
    )
    rpsbml_json = os_path.join(
        Main_rplibs.data_path,
        'rpsbml.json'
    )
    ref_name  = 'RetroPath_Pathway_1_1'
    ref_score = 0.5684564101634014


    def setUp(self):
        super().setUp()
        # objects below have to be created for each test instance
        # since some tests can modified them
        self.rpsbml  = rpSBML(
            inFile = self.rpsbml_path,
            logger = self.logger
        )
        self.gem = rpSBML(
            inFile = extract_gz(
                self.gem_path_gz,
                self.temp_d
            )
        )


    def test_initEmpty(self):
        rpSBML(name='rpSBML_test', logger=self.logger)


    def test_initWithInFile(self):
        self.assertEqual(self.rpsbml.getName(), self.ref_name)


    def test_initWithDocument(self):
        rpsbml = rpSBML(rpsbml = self.rpsbml)
        self.assertEqual(rpsbml.getName(), self.ref_name)


    def test_initWithModelName(self):
        rpsbml = rpSBML(name=self.rpsbml.getName())
        self.assertEqual(rpsbml.getName(), self.ref_name)


    def test_initWithNothing(self):
        rpsbml = rpSBML()
        self.assertEqual(rpsbml.getName(), 'dummy')


    def test_score(self):
        self.rpsbml.compute_score()
        self.assertAlmostEqual(self.rpsbml.getScore(), self.ref_score)


    def test_globalscore(self):
        global_score = self.rpsbml.compute_globalscore()
        self.assertAlmostEqual(global_score, 0.5760957019721074)


    '''
    def test_dictRPpathway(self):
        self.assertDictEqual(self.rpsbml._dictRPpathway(), self.data['dictrppathway'])
    '''


    def test_nameToSbmlId(self):
        self.assertEqual(
            self.rpsbml._nameToSbmlId('test123-_!"£$%^&*(){}@~><>?'),
            'test123___________________'
        )


    def test_genMetaID(self):
        self.assertEqual(
            self.rpsbml._genMetaID('test123'),
            'ecd71870d1963316a97e3ac3408c9835ad8cf0f3c1bc703527c30265534f75ae'
        )


    def test_toDict(self):
        with open(self.rpsbml_json, 'r') as f:
            self.assertDictEqual(
                self.rpsbml.toDict(),
                json_load(f)
            )


    def test_toJSON(self):
        with open(self.rpsbml_json, 'r') as f:
            self.assertDictEqual(
                json_loads(self.rpsbml.toJSON()),
                json_load(f)
            )


    def test_updateBRSynthPathway(self):
        rpsbml = rpSBML(
            inFile = os_path.join(
                self.data_path,
                'rpsbml_empty_pathway.xml'
            )
        )
        self_rpsbml_dict = self.rpsbml.toDict()
        rpsbml.updateBRSynthPathway(self_rpsbml_dict)
        self.assertDictEqual(
            rpsbml.toDict(),
            self_rpsbml_dict
        )


    def test_readRPrules(self):
        infile = os_path.join(
            self.data_path,
            'rp_rules.json'
        )
        with open(infile, 'r') as f:
            self.assertDictEqual(
                self.rpsbml.readRPrules(),
                json_loads(f.read())
            )


    def test_readRPspecies(self):
        infile = os_path.join(
            self.data_path,
            'rp_species.json'
        )
        with open(infile, 'r') as f:
            self.assertDictEqual(
                self.rpsbml.readRPspecies(),
                json_loads(f.read())
            )


    def test_readRPpathwayIDs(self):
        self.assertCountEqual(
            self.rpsbml.readGroupMembers('rp_pathway'),
            ['rxn_1', 'rxn_2', 'rxn_3']
        )


    def test_readUniqueRPspecies(self):
        self.assertCountEqual(
            self.rpsbml.readUniqueRPspecies(),
            [
                'TARGET_0000000001__64__MNXC3',
                'MNXM13__64__MNXC3',
                'CMPD_0000000004__64__MNXC3',
                'MNXM1__64__MNXC3',
                'MNXM20__64__MNXC3',
                'CMPD_0000000013__64__MNXC3',
                'MNXM89557__64__MNXC3',
                'MNXM5__64__MNXC3',
                'MNXM7__64__MNXC3',
                'MNXM9__64__MNXC3',
                'MNXM6__64__MNXC3',
                'MNXM3__64__MNXC3'
            ]
        )


    def test_speciesExists(self):
        self.assertTrue(self.rpsbml.speciesExists('MNXM89557'))
        self.assertFalse(self.rpsbml.speciesExists('test'))


    def test_isSpeciesProduct(self):
        self.assertTrue(self.rpsbml.isSpeciesProduct('TARGET_0000000001__64__MNXC3'))
        self.assertFalse(self.rpsbml.isSpeciesProduct('MNXM1__64__MNXC3'))


    def test_compareBRSYNTHAnnotations(self):
        # self.assertTrue(self.rpsbml.compareBRSYNTHAnnotations(
        #     self.rpsbml.getModel().getSpecies('MNXM89557__64__MNXC3').getAnnotation(),
        #     self.rpsbml.getModel().getSpecies('MNXM89557__64__MNXC3').getAnnotation()))
        # self.assertFalse(self.rpsbml.compareBRSYNTHAnnotations(
        #     self.rpsbml.getModel().getSpecies('MNXM89557__64__MNXC3').getAnnotation(),
        #     self.rpsbml.getModel().getSpecies('CMPD_0000000013__64__MNXC3').getAnnotation()))
        self.assertFalse(self.gem.compareBRSYNTHAnnotations(
            self.gem.getModel().getSpecies('M_2pg_c').getAnnotation(),
            self.gem.getModel().getSpecies('M_13dpg_c').getAnnotation()))


    def test_compareMIRIAMAnnotations(self):
        self.assertTrue(self.rpsbml.compareMIRIAMAnnotations(
            self.rpsbml.getModel().getSpecies('MNXM89557__64__MNXC3').getAnnotation(),
            self.rpsbml.getModel().getSpecies('MNXM89557__64__MNXC3').getAnnotation()))
        self.assertFalse(self.rpsbml.compareMIRIAMAnnotations(
            self.rpsbml.getModel().getSpecies('MNXM89557__64__MNXC3').getAnnotation(),
            self.rpsbml.getModel().getSpecies('CMPD_0000000013__64__MNXC3').getAnnotation()))


    def test_createReturnFluxParameter(self):
        #return feature
        param = self.rpsbml.createReturnFluxParameter(
            None,
            parameter_id = 'B_999999'
        )
        self.assertEqual(param.id, 'B_999999')
        self.assertEqual(param.value, 999999.0)
        #create feature
        new = rpSBML(name='test', logger=self.logger)
        new.createModel('test_name', 'test_id')
        param = new.createReturnFluxParameter(8888.0)
        self.assertEqual(param.id, 'B_8888_0')
        self.assertEqual(param.value, 8888.0)


    def test_readMIRIAMAnnotation(self):
        self.assertDictEqual(self.rpsbml.readMIRIAMAnnotation(
            self.rpsbml.getModel().getReaction('rxn_3').getAnnotation()),
            {'ec-code': ['4.1.1.17']}
        )
        self.assertDictEqual(self.gem.readMIRIAMAnnotation(
            self.gem.getModel().getReaction('R_ALATA_D2').getAnnotation()),
            {
                'bigg':     ['ALATA_D2'],
                'biocyc':   ['RXN0-5240'],
                'kegg':     ['R01147'],
                'metanetx': ['MNXR95697'],
                'rhea':     ['28562', '28563', '28564', '28565']
            }
        )


    def test_readReactionSpecies(self):
        self.assertDictEqual(self.rpsbml.readReactionSpecies(
            self.rpsbml.getModel().getReaction('rxn_3')),
                {
                    'left': {
                        'CMPD_0000000004__64__MNXC3': 1,
                        'MNXM1__64__MNXC3': 1
                    },
                    'right': {
                        'TARGET_0000000001__64__MNXC3': 1,
                        'MNXM13__64__MNXC3': 1
                    }
                }
        )


    def test_mergeModels(self):
        target_rpsbml = rpSBML(
            inFile = self.e_coli_model_path,
            logger = self.logger
        )
        merged_rpsbml, reactions_in_both = rpSBML.mergeModels(
            self.rpsbml,
            target_rpsbml,
            logger = self.logger
            )
        ref_target_rpsbml = rpSBML(
            inFile = self.merged_path,
            logger = self.logger
        )
        self.assertEqual(merged_rpsbml, ref_target_rpsbml)


    def test_mergeFiles(self):
        with NamedTemporaryFile() as tempf:
            rpSBML.mergeFiles(
                self.rpsbml_path,
                self.e_coli_model_path,
                tempf.name,
                logger = self.logger
            )
            self.assertEqual(
                rpSBML(
                    inFile = tempf.name,
                    logger = self.logger
                ),
                rpSBML(
                    inFile = self.merged_path,
                    logger = self.logger
                )
            )


    def test_print(self):
        comp_xref_path = os_path.join(
            Main_rplibs.data_path,
            'comp_xref.json'
        )
        # comp_xref = {'bigg': ['c', 'c_c'], 'mnx': ['MNXC3'], 'cco': ['cco-clath-end-ves', 'cco-clath-ves', 'cco-coated-ves', 'cco-copii-ves', 'cco-cytoplasm', 'cco-cytosol', 'cco-early-end-lum', 'cco-early-endo', 'cco-early-endo-mem', 'cco-endo-lum', 'cco-endo-mem', 'cco-endocyt-ves', 'cco-endosome', 'cco-er-golgi-ves', 'cco-late-end-lum', 'cco-late-endo', 'cco-late-endo-mem', 'cco-mbody-mem', 'cco-micro-lum', 'cco-microbody', 'cco-ribosome', 'cco-sec-granule', 'c'], 'go': ['0000015', '0000133', '0000145', '0000153', '0000164', '0000177', '0000229', '0000235', '0000242', '0000308', '0000345', '0000407', '0000787', '0000789', '0000809', '0000924', '0000927', '0000930', '0000931', '0000932', '0002079', '0002081', '0002096', '0002144', '0002186', '0005737', '0005768', '0005769', '0005770', '0005771', '0005786', '0005793', '0005814', '0005821', '0005822', '0005823', '0005824', '0005825', '0005829', '0005831', '0005832', '0005833', '0005835', '0005836', '0005840', '0005850', '0005851', '0005852', '0005853', '0005854', '0005859', '0005863', '0005881', '0005883', '0005938', '0005945', '0005946', '0005948', '0005950', '0005951', '0005960', '0005964', '0005965', '0005968', '0005969', '0005971', '0008074', '0008274', '0008275', '0008303', '0008352', '0008385', '0009316', '0009317', '0009320', '0009321', '0009328', '0009329', '0009331', '0009332', '0009333', '0009336', '0009339', '0009340', '0009343', '0009345', '0009346', '0009348', '0009357', '0009359', '0009361', '0009376', '0009382', '0009504', '0009524', '0009525', '0009574', '0010005', '0010009', '0010169', '0010494', '0012507', '0014705', '0015934', '0015935', '0016281', '0016282', '0016465', '0016528', '0017102', '0017109', '0018444', '0019037', '0019197', '0019812', '0019813', '0020026', '0022625', '0022626', '0022627', '0030016', '0030017', '0030018', '0030077', '0030078', '0030079', '0030080', '0030081', '0030082', '0030117', '0030118', '0030119', '0030120', '0030123', '0030124', '0030125', '0030127', '0030128', '0030131', '0030133', '0030134', '0030135', '0030136', '0030139', '0030141', '0030142', '0030143', '0030314', '0030485', '0030486', '0030658', '0030659', '0030661', '0030662', '0030665', '0030666', '0030667', '0030668', '0030669', '0030670', '0030671', '0030700', '0030877', '0030891', '0030929', '0030930', '0030931', '0031045', '0031082', '0031083', '0031084', '0031085', '0031088', '0031089', '0031091', '0031092', '0031093', '0031201', '0031209', '0031250', '0031251', '0031313', '0031410', '0031414', '0031415', '0031416', '0031417', '0031430', '0031500', '0031515', '0031560', '0031561', '0031562', '0031592', '0031597', '0031600', '0031603', '0031606', '0031609', '0031612', '0031615', '0031633', '0031672', '0031673', '0031674', '0031680', '0031901', '0031902', '0031903', '0031904', '0031905', '0031906', '0031907', '0032009', '0032019', '0032047', '0032068', '0032117', '0032127', '0032449', '0032477', '0032578', '0032585', '0032593', '0032838', '0032839', '0032921', '0033093', '0033095', '0033107', '0033110', '0033116', '0033117', '0033118', '0033150', '0033162', '0033254', '0033282', '0033290', '0033291', '0033391', '0033596', '0033675', '0034045', '0034081', '0034082', '0034083', '0034270', '0034274', '0034422', '0034430', '0034451', '0034466', '0034467', '0034492', '0034493', '0034515', '0034519', '0034709', '0034715', '0034719', '0034730', '0034731', '0034741', '0034743', '0034744', '0034745', '0034746', '0034748', '0034749', '0034750', '0034752', '0034774', '0034777', '0034973', '0034996', '0035550', '0035579', '0035580', '0036000', '0036007', '0036186', '0036379', '0036396', '0036411', '0036457', '0036464', '0036501', '0042151', '0042470', '0042566', '0042579', '0042581', '0042583', '0042584', '0042587', '0042588', '0042589', '0042599', '0042709', '0042716', '0042717', '0042718', '0042735', '0042788', '0042827', '0043034', '0043186', '0043223', '0043265', '0043291', '0043292', '0043293', '0043527', '0043528', '0043540', '0043597', '0043598', '0043600', '0043614', '0043659', '0043660', '0043661', '0043662', '0043698', '0043699', '0043700', '0043701', '0043702', '0044207', '0044222', '0044223', '0044227', '0044310', '0044312', '0044352', '0044353', '0044354', '0044391', '0044433', '0044438', '0044440', '0044444', '0044445', '0044448', '0044449', '0044450', '0044599', '0044797', '0044840', '0044841', '0044842', '0045009', '0045169', '0045170', '0045239', '0045243', '0045246', '0045247', '0045248', '0045249', '0045250', '0045251', '0045253', '0045254', '0045293', '0045334', '0045335', '0045336', '0045495', '0046859', '0048471', '0048492', '0048494', '0048500', '0048501', '0048770', '0055028', '0055037', '0055038', '0055087', '0055120', '0060198', '0060199', '0060200', '0060201', '0060202', '0060203', '0060204', '0060205', '0060293', '0060417', '0060418', '0060473', '0061200', '0061201', '0061202', '0061493', '0061496', '0061497', '0061498', '0061499', '0061645', '0061673', '0061702', '0061802', '0061803', '0070032', '0070033', '0070044', '0070045', '0070046', '0070047', '0070048', '0070049', '0070065', '0070066', '0070067', '0070068', '0070081', '0070082', '0070083', '0070088', '0070111', '0070112', '0070113', '0070114', '0070115', '0070116', '0070117', '0070118', '0070214', '0070319', '0070355', '0070356', '0070381', '0070382', '0070422', '0070441', '0070554', '0070695', '0070766', '0070768', '0070820', '0070821', '0070826', '0070937', '0070992', '0070993', '0071075', '0071203', '0071212', '0071254', '0071341', '0071439', '0071513', '0071521', '0071540', '0071541', '0071546', '0071547', '0071563', '0071598', '0071682', '0071818', '0072379', '0072380', '0072516', '0072557', '0072558', '0072559', '0080085', '0090569', '0090651', '0090652', '0090653', '0090654', '0090726', '0090730', '0097013', '0097169', '0097208', '0097209', '0097232', '0097233', '0097234', '0097342', '0097361', '0097422', '0097433', '0097443', '0097486', '0097487', '0097488', '0097489', '0097495', '0097512', '0097547', '0097632', '0097633', '0097634', '0097654', '0097721', '0098538', '0098539', '0098559', '0098560', '0098566', '0098593', '0098594', '0098723', '0098832', '0098898', '0098922', '0098956', '0099073', '0099078', '0099501', '0099503', '0099522', '0099523', '0099524', '0099568', '0099738', '0101002', '0101003', '0106002', '0120002', '1902502', '1902507', '1903144', '1903145', '1903754', '1904115', '1904511', '1904564', '1904724', '1904813', '1905720', '1905721', '1990005', '1990006', '1990008', '1990061', '1990070', '1990133', '1990140', '1990198', '1990200', '1990201', '1990202', '1990205', '1990220', '1990221', '1990228', '1990229', '1990230', '1990231', '1990232', '1990233', '1990235', '1990252', '1990257', '1990316', '1990329', '1990330', '1990334', '1990457', '1990463', '1990499', '1990500', '1990565', '1990574', '1990588', '1990633', '1990657', '1990658', '1990673', '1990701', '1990723', '1990728', '1990730', '1990733', '1990753', '1990783', '1990788', '1990917'], 'seed': ['c', 'c0', 'cytosol'], 'name': ['cytosol']}
        with open(comp_xref_path, 'r') as f:
            comp_xref = json_load(f)
        rpsbml = rpSBML(
            name='rpSBML_test',
            logger=self.logger
        )
        rpsbml.genericModel(
            'RetroPath_Pathway_test',
            'RP_model_test',
            comp_xref,
            'MNXC3',
            999999,
            0
        )
        rpsbml.createGroup('rp_pathway')
        rpsbml.createGroup('central_species')
        with NamedTemporaryFile() as tempf:
            rpsbml.writeSBML(tempf.name)
            self.assertListEqual(
                list(
                    io_open(tempf.name)
                ),
                list(
                    io_open(
                        os_path.join(
                            self.data_path,
                            'rpSBML_test_sbml.xml'
                        )
                    )
                )
            )
