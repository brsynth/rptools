from unittest import TestCase
from os import path as os_path
from tempfile import (
    TemporaryDirectory,
    mkdtemp
)
from shutil import rmtree
from cobra.core.solution import Solution  as cobra_solution
from json import load as json_load
from brs_utils import (
    create_logger,
    extract_gz
)
from rptools.rpfba.rpFBA import (
    runFBA
)
from rptools.rplibs import rpSBML


class Test_rpFBA(TestCase):


    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    e_coli_model_path_gz = os_path.join(
        data_path,
        'e_coli_model.sbml.gz'
    )
    pathway_path = os_path.join(
        data_path,
        'pathway.json'
    )

    def test(self):
        with open(self.pathway_path, 'r') as fp:
            pathway = json_load(fp)
        results = runFBA(
            pathway=pathway,
            gem_sbml_path=self.e_coli_model_path_gz,
            sim_type='fraction'
        )
        self.assertDictEqual(
            results,
            {
                "reactions": {
                    "rxn_1": 1.3648925522849882,
                    "rxn_2": 1.3648925522849882,
                    "rxn_3": 1.3648925522849882,
                    "rxn_4": 1.3648925522849882
                },
                "pathway": {
                    "biomass": 0.57290585662576
                }
            }
        )


    # merged_path_gz = os_path.join(
    #     data_path,
    #     'merged.xml.gz'
    # )
    # nb_rpsbml_paths = 2
    # # rxn_tgt = 'rxn_target'
    # # pathway_id = 'rp_pathway'
    # fraction_scores = [
    #     [
    #         2.3076923076923888,
    #         3.6794124272706443,
    #         3.6794124272706443
    #     ],
    #     [
    #         1.3296695186776557,
    #         0.7638744755010182,
    #         0.7638744755010182
    #     ]
    # ]
    # fba_scores = [
    #         9.230769230769237,
    #         5.144315545243618,
    #         3.398422535211268
    #     ]
    # pfba_scores = [
    #         859.3846153846168,
    #         471.6390839533362,
    #         382.1106535211269
    #     ]

    # def setUp(self):
    #     self.logger = create_logger(__name__, 'ERROR')

    #     self.rpsbml_paths = [
    #         os_path.join(
    #             self.data_path,
    #             'rpsbml_'+str(i+1)+'.xml'
    #         ) for i in range(self.nb_rpsbml_paths)
    #     ]

    #     self.rpsbmls = [
    #         rpSBML(
    #             inFile = self.rpsbml_paths[i],
    #             logger = self.logger
    #         ) for i in range(self.nb_rpsbml_paths)
    #     ]

    #     with TemporaryDirectory() as temp_d:
    #         self.e_coli_model_path = extract_gz(
    #             self.e_coli_model_path_gz,
    #             temp_d
    #         )
    #         self.merged_path = extract_gz(
    #             self.merged_path_gz,
    #             temp_d
    #         )
    #         self.rpsbml_model = rpSBML(
    #             inFile = self.e_coli_model_path,
    #             logger = self.logger
    #         )
    #         self.merged_rpsbml_1 = rpSBML(
    #             inFile = self.merged_path,
    #             logger = self.logger
    #         )
    #         self.merged_rpsbml_2, reactions_in_both = rpSBML.mergeModels(
    #             source_rpsbml = self.rpsbmls[0],
    #             target_rpsbml = self.rpsbml_model,
    #             logger = self.logger
    #         )
    #         self.merged_rpsbml_3, reactions_in_both = rpSBML.mergeModels(
    #             source_rpsbml = self.rpsbmls[1],
    #             target_rpsbml = self.rpsbml_model,
    #             logger = self.logger
    #         )

    # def test_fba(self):
    #     objective_id = 'obj_' + self.rxn_tgt
    #     for i in range(len(self.fba_scores)):
    #         with self.subTest(
    #             i = i
    #         ):
    #             ref_score = self.fba_scores[i]
    #             rpsbml = getattr(
    #                 self,
    #                 'merged_rpsbml_' + str(i+1)
    #             )

    #             objective_id = rpsbml.find_or_create_objective(
    #                 reactions = [self.rxn_tgt],
    #                 coefficients = [1.0],
    #                 is_max = True,
    #                 objective_id = objective_id
    #             )

    #             rpsbml.search_isolated_species()

    #             cobra_solution = rp_fba(
    #                     rpsbml = rpsbml,
    #                     objective_id = objective_id,
    #                     logger = self.logger
    #             )

    #             self._test(
    #                 ref_score,
    #                 rpsbml,
    #                 cobra_solution,
    #                 objective_id
    #             )

    # def test_pfba(self):
    #     objective_id = 'obj_' + self.rxn_tgt
    #     for i in range(len(self.pfba_scores)):
    #         with self.subTest(
    #             i = i
    #         ):
    #             ref_score = self.pfba_scores[i]
    #             rpsbml = getattr(
    #                 self,
    #                 'merged_rpsbml_' + str(i+1)
    #             )

    #             objective_id = rpsbml.find_or_create_objective(
    #                 reactions = [self.rxn_tgt],
    #                 coefficients = [1.0],
    #                 is_max = True,
    #                 objective_id = objective_id
    #             )

    #             rpsbml.search_isolated_species()

    #             cobra_solution = rp_pfba(
    #                     rpsbml = rpsbml,
    #                     objective_id = objective_id,
    #                     logger = self.logger
    #             )

    #             self._test(
    #                 ref_score,
    #                 rpsbml,
    #                 cobra_solution,
    #                 objective_id
    #             )

    # def _test(
    #     self,
    #     ref_score: float,
    #     rpsbml: rpSBML,
    #     cobra_solution: cobra_solution,
    #     objective_id: str
    # ) -> None:
    #     self.assertTrue(rpsbml)
    #     self.assertAlmostEqual(
    #         cobra_solution.objective_value,
    #         ref_score
    #     )

    #     write_results(
    #         rpsbml = rpsbml,
    #         objective_id = objective_id,
    #         cobra_results = cobra_solution,
    #         pathway_id = self.pathway_id,
    #         logger = self.logger
    #     )
    #     # make sure that the results are written to the file
    #     pathway = rpsbml.toDict()['pathway']['brsynth']
    #     self.assertAlmostEqual(
    #         pathway['fba_obj_'+self.rxn_tgt]['value'],
    #         ref_score
    #     )

    # def test_fraction(self):
    #     for i in range(len(self.fraction_scores)):
    #         scores = self.fraction_scores[i]
    #         rpsbml = getattr(
    #             self,
    #             'merged_rpsbml_' + str(i+1)
    #         )
    #         rpsbml.search_isolated_species()
    #         cobra_results, rpsbml = rp_fraction(
    #                 rpsbml = rpsbml,
    #             src_rxn_id = 'biomass',
    #             src_coeff = 1.0,
    #             tgt_rxn_id = self.rxn_tgt,
    #             tgt_coeff = 1.0,
    #             frac_of_src = 0.75,
    #                 is_max = True,
    #             pathway_id = self.pathway_id,
    #             objective_id = None,
    #                 logger = self.logger
    #         )

    #         self.assertTrue(rpsbml)
    #         self.assertAlmostEqual(
    #             cobra_results.objective_value,
    #             scores[0]
    #         )

    #         # make sure that the results are written to the file
    #         pathway = rpsbml.toDict()['pathway']['brsynth']
    #         self.assertAlmostEqual(
    #             pathway['fba_obj_'+self.rxn_tgt+'__restricted_biomass']['value'],
    #             scores[0]
    #         )
    #         self.assertAlmostEqual(
    #             pathway['fba_obj_biomass']['value'],
    #             scores[1]
    #         )
