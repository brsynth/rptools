from glob import glob
from zipfile import ZipFile
from types import SimpleNamespace
from tempfile import NamedTemporaryFile

from os import path as os_path

from rptools.rpfba.rpfba import preprocess, runFBA
from rptools.rplibs import rpPathway, rpSBML
from main_rpfba import Main_rpfba


class Test_rpFBA(Main_rpfba):
    def setUp(self):
        super().setUp()
        input_zip = ZipFile(self.cr_path)
        input_zip.extractall(path=os_path.join(self.temp_d, "cr_fba"))
        output_zip = ZipFile(self.fba_path)
        output_zip.extractall(path=os_path.join(self.temp_d, "lycopene_fba"))

    def tearDown(self):
        super().tearDown()

    def test_runFBA(self):
        # TODO: taking accound extra args, like medium
        def _extract_var(dirname):
            files = glob(os_path.join(dirname, "*xml"))
            basenames = [os_path.basename(x) for x in files]
            names = [x.split(".")[0] for x in basenames]
            sims = [x.split(".")[1] for x in basenames]
            assert len(files) == len(names) == len(sims)
            return (files, names, sims)

        def _extract_res_from_file(filename):
            rp_pathway = rpPathway(infile=filename)
            res = {}
            res["pathway"] = rp_pathway.get_fba()
            res["reactions"] = {}
            for rid in rp_pathway.get_reactions():
                res["reactions"][rid] = rp_pathway.get_reaction(rid).get_fba()
            return res

        def _format_dict(old, new={}):
            print("old", old, "new", new)
            for k, v in old.items():
                if isinstance(v, dict):
                    new[k] = _format_dict(old.get(k, {}), v)
                else:
                    if isinstance(v, str):
                        new[k] = v
                    else:
                        new[k] = round(float(v), 2)
            return new

        files, names, sims = _extract_var(os_path.join(self.temp_d, "lycopene_fba"))

        args = SimpleNamespace(
            model_file=self.e_coli_model_path,
            compartment_id="c",
            biomass_rxn_id="biomass",
            objective_rxn_id="rxn_target",
            ignore_orphan_species=True,
            sim="fraction",
            fraction_of=0.75,
            merge=True
        )
        for ix in range(len(files)):
            args.pathway_file = os_path.join(self.temp_d, "cr_fba", names[ix] + ".xml")
            (
                merged_model,
                pathway,
                ids
            ) = preprocess(args=args)
            results = runFBA(
                model=merged_model,
                compartment_id=ids['comp_id'],
                biomass_rxn_id=ids['biomass_rxn_id'],
                objective_rxn_id=ids['obj_rxn_id'],
                sim_type=sims[ix],
                fraction_coeff=args.fraction_of
            )

            res_previous = _format_dict(_extract_res_from_file(files[ix]))
            res_run_fba = _format_dict(
                {x: y for x, y in results.items() if x in ["pathway", "reactions"]}
            )

            self.assertDictEqual(res_previous, res_run_fba)
