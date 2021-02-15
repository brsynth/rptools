from unittest import TestCase
from os import path as os_path
from rptools.rpfba import rpFBA
from rptools.rpfba.rpFBA import (
    # rp_fba,
    rp_fraction,
    # rp_pfba,
    write_results
)
from rptools.rplibs import rpSBML
from tempfile import (
    TemporaryDirectory,
    mkdtemp
)
from shutil    import rmtree
from brs_utils import (
    create_logger,
    extract_gz
)


class Test_rpFBA(TestCase):


    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    merged_path_gz = os_path.join(
        data_path,
        'merged.xml.gz'
    )
    rpsbml_path = os_path.join(
        data_path,
        'rp_1_1_sbml.xml'
    )
    e_coli_model_path_gz = os_path.join(
        data_path,
        'e_coli_model.sbml.gz'
    )

    rxn_tgt = 'rxn_target'
    pathway_id = 'rp_pathway'


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

        self.rpsbml = rpSBML(
            inFile = self.rpsbml_path,
            logger = self.logger
        )

        with TemporaryDirectory() as temp_d:
            self.e_coli_model_path = extract_gz(
                self.e_coli_model_path_gz,
                temp_d
            )
            self.merged_path = extract_gz(
                self.merged_path_gz,
                temp_d
            )
            self.rpsbml_model = rpSBML(
                inFile = self.e_coli_model_path,
                logger = self.logger
            )
            self.merged_rpsbml_1 = rpSBML(
                inFile = self.merged_path,
                logger = self.logger
            )
            self.merged_rpsbml_2, reactions_in_both = rpSBML.mergeModels(
                source_rpsbml = self.rpsbml,
                target_rpsbml = self.rpsbml_model,
                logger = self.logger
            )


    def test_fba_pfba(self):
        objective_id = 'obj_' + self.rxn_tgt
        ref_scores = {
            'fba': [
                9.230769230769237,
                5.144315545243618
            ],
            'pfba': [
                859.3846153846168,
                471.6390839533362
            ]
        }
        for f, scores in ref_scores.items():
            func = getattr(rpFBA, 'rp_' + f)
            for i in range(len(scores)):
                with self.subTest(
                    func = func,
                    scores = scores,
                    i = i
                ):
                    ref_score = scores[i]
                    rpsbml = getattr(
                        self,
                        'merged_rpsbml_' + str(i+1)
                    )

                    objective_id = rpsbml.find_or_create_objective(
                        reactions = [self.rxn_tgt],
                        coefficients = [1.0],
                        is_max = True,
                        objective_id = objective_id
                    )
                    rpsbml.activateObjective(
                        objective_id = objective_id,
                        plugin = 'fbc'
                    )

                    cobra_solution = func(
                            rpsbml = rpsbml,
                            logger = self.logger
                    )

                    self.assertTrue(rpsbml)
                    self.assertAlmostEqual(
                        cobra_solution.objective_value,
                        ref_score
                    )

                    write_results(
                        rpsbml = rpsbml,
                        objective_id = objective_id,
                        cobra_results = cobra_solution,
                        pathway_id = self.pathway_id,
                        logger = self.logger
                    )
                    # make sure that the results are written to the file
                    pathway = rpsbml.toDict()['pathway']['brsynth']
                    self.assertAlmostEqual(
                        pathway['fba_obj_'+self.rxn_tgt]['value'],
                        ref_score
                    )


    def test_fraction(self):
        ref_scores = [
            [
                2.3076923076923888,
                3.6794124272706443
            ],
            [
                1.3296695186776557,
                0.7638744755010182
            ]
        ]

        for i in range(len(ref_scores)):
            scores = ref_scores[i]
            rpsbml = getattr(
                self,
                'merged_rpsbml_' + str(i+1)
            )
            cobra_results, rpsbml = rp_fraction(
                    rpsbml = rpsbml,
                src_rxn_id = 'biomass',
                src_coeff = 1.0,
                tgt_rxn_id = self.rxn_tgt,
                tgt_coeff = 1.0,
                frac_of_src = 0.75,
                    is_max = True,
                pathway_id = self.pathway_id,
                objective_id = None,
                    logger = self.logger
            )

            self.assertTrue(rpsbml)
            self.assertAlmostEqual(
                cobra_results.objective_value,
                scores[0]
            )

            # make sure that the results are written to the file
            pathway = rpsbml.toDict()['pathway']['brsynth']
            self.assertAlmostEqual(
                pathway['fba_obj_'+self.rxn_tgt+'__restricted_biomass']['value'],
                scores[0]
            )
            self.assertAlmostEqual(
                pathway['fba_obj_biomass']['value'],
                scores[1]
            )
