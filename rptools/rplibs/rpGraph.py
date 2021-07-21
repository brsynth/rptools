import networkx as nx
from networkx.readwrite import json_graph
import logging
import os
import itertools
import numpy as np
import random
from typing import(
    List,
    Callable
)


class rpGraph:
    """The class that hosts the networkx related functions
    """
    def __init__(self,
                 rpsbml,
                 is_gem_sbml,
                 pathway_id,
                 central_species_group_id,
                 sink_species_group_id,
                 logger=logging.getLogger(__name__)):
        """Constructor of the class

        Automatically constructs the network when calling the construtor

        :param rpsbml: The rpSBML object
        :param pathway_id: The pathway id of the heterologous pathway
        :param species_group_id: The id of the central species

        :type rpsbml: rpSBML
        :type pathway_id: str
        :type species_group_id: str
        """
        self.logger = logger

        self.logger.debug('New instance of rpGraph')

        self.rpsbml = rpsbml
        self.pathway_id = pathway_id
        self.central_species_group_id = central_species_group_id
        self.sink_species_group_id = sink_species_group_id
        self.G = None
        self.pathway_id = pathway_id
        self.num_reactions = 0
        self.num_species = 0
        if rpsbml:
            self.__makeGraph(is_gem_sbml, pathway_id, central_species_group_id, sink_species_group_id)


    ######################################################################################################
    ######################################### Private Function ###########################################
    ######################################################################################################


    ################################# Analyse and make graph #####################


    #TODO: add the compartments to the species and reactions node descriptions
    def __makeGraph(
        self,
        is_gem_sbml,
        pathway_id,
        central_species_group_id,
        sink_species_group_id
    ):
        """Private function that constructs the networkx graph

        :param is_gem_sbml: Determine what type of graph to build. If True then all the species and reactions will be added and not just the heterologous pathway.
        :param pathway_id: The pathway id of the heterologous pathway
        :param species_group_id: The id of the central species

        :type pathway_id: str
        :type species_group_id: str

        :return: None
        :rtype: None
        """

        rpsbml_model = self.rpsbml.getModel()
        # rp_species = [rpsbml_model.getSpecies(i) for i in self.rpsbml.readUniqueRPspecies(pathway_id)]
        groups = rpsbml_model.getPlugin('groups')
        # c_s = self.rpsbml.getGroup(central_species_group_id)
        s_s = self.rpsbml.getGroup(sink_species_group_id)
        # rp_central_species_id = [i.getIdRef() for i in c_s.getListOfMembers()]
        rp_sink_species_id = [i.getIdRef() for i in s_s.getListOfMembers()]
        rp_pathway = self.rpsbml.getGroup(pathway_id)
        rp_species_id = self.rpsbml.readUniqueRPspecies()
        rp_reactions_id = [i.getIdRef() for i in rp_pathway.getListOfMembers()]
        self.logger.debug('rp_reactions_id: '+str(rp_reactions_id))
        self.G = nx.DiGraph(
            brsynth = self.rpsbml.readBRSYNTHAnnotation(
                annot = rp_pathway.getAnnotation(),
                logger = self.logger
            )
        )

        #### add ALL the species and reactions ####
        # nodes
        for species in rpsbml_model.getListOfSpecies():
            is_central = False
            is_sink = False
            is_rp_pathway = False
            if species.getId() in rp_species_id:
                is_rp_pathway = True
            # if species.getId() in rp_central_species_id:
            #     is_central = True
            if species.getId() in rp_sink_species_id:
                is_sink = True
            # add it if GEM then all, or if rp_pathway
            if is_rp_pathway or is_gem_sbml:
                self.num_species += 1
                self.G.add_node(
                    species.getId(),
                    type = 'species',
                    name = species.getName(),
                    miriam = self.rpsbml.readMIRIAMAnnotation(
                        species.getAnnotation()
                    ),
                    brsynth = self.rpsbml.readBRSYNTHAnnotation(
                        annot = species.getAnnotation(),
                        logger = self.logger
                    ),
                    rp_trunk_species = is_central,
                    rp_sink_species = is_sink,
                    rp_pathway = is_rp_pathway
                )

        for reaction in rpsbml_model.getListOfReactions():
            is_rp_pathway = False
            if reaction.getId() in rp_reactions_id:
                is_rp_pathway = True
            if is_rp_pathway or is_gem_sbml:
                self.num_reactions += 1
                self.G.add_node(
                    reaction.getId(),
                    type = 'reaction',
                    miriam = self.rpsbml.readMIRIAMAnnotation(
                        reaction.getAnnotation()
                    ),
                    brsynth = self.rpsbml.readBRSYNTHAnnotation(
                        reaction.getAnnotation(),
                        logger = self.logger
                    ),
                    rp_pathway = is_rp_pathway
                )

        # edges
        for reaction in rpsbml_model.getListOfReactions():
            self.logger.debug('Adding edges for the reaction: '+str(reaction.getId()))
            if reaction.getId() in rp_reactions_id or is_gem_sbml:
                for reac in reaction.getListOfReactants():
                    # if reac.species == 'TARGET_0000000001':
                    #     self.logger.info('\taAdding edge '+str(reac.species)+' --> '+str(reaction.getId()))
                    self.logger.debug('\taAdding edge '+str(reac.species)+' --> '+str(reaction.getId()))
                    self.G.add_edge(
                        reac.species,
                        reaction.getId(),
                        stoichio = reac.stoichiometry
                    )
                for prod in reaction.getListOfProducts():
                    self.logger.debug('\taAdding edge '+str(reaction.getId())+' --> '+str(prod.species))
                    self.G.add_edge(
                        reaction.getId(),
                        prod.species,
                        stoichio=reac.stoichiometry
                    )


    def __isolatedSpecies(
        self,
        func_eq_0: Callable,
        func_gt_0: Callable,
        species: List[str] = [],
        only_central: bool = False,
        only_rp_pathway: bool = True
    ) -> List[str]:
        """Private function that returns the single parent species that are consumed only

        :param only_central: Focus on the central species only

        :type only_central: bool

        :return: List of node ids
        :rtype: list
        """
        _species = []
        if not species:
            species = self.G.nodes()
        for spe_id in species:
        # for node_name in self.G.nodes():
            node = self.G.nodes.get(spe_id)
            if node['type'] == 'species':
                if node['rp_pathway']:
                    # NOTE: if central species then must also be rp_pathway species
                    if (
                        (
                            only_central
                            and node['rp_trunk_species']
                        )
                        or (
                            only_rp_pathway
                            and node['rp_pathway']
                        )
                        or (
                            not only_central
                            and not only_rp_pathway
                        )
                    ):
                        self.logger.debug(spe_id)
                        self.logger.debug(f'succ: {list(func_gt_0(spe_id))}')
                        self.logger.debug(f'pred: {list(func_eq_0(spe_id))}')
                        if (
                            len(list(func_gt_0(spe_id))) > 0
                            and len(list(func_eq_0(spe_id))) == 0
                        ):
                            _species.append(spe_id)

        return _species


    def onlyConsumedSpecies(self, species=[], only_central=False, only_rp_pathway=True):
        """Private function that returns the single parent species that are consumed only

        :param only_central: Focus on the central species only

        :type only_central: bool

        :return: List of node ids
        :rtype: list
        """
        return self.__isolatedSpecies(
            species=species,
            only_central=only_central,
            only_rp_pathway=only_rp_pathway,
            func_eq_0=self.G.predecessors,
            func_gt_0=self.G.successors
        )


    def onlyProducedSpecies(self, species=[], only_central=False, only_rp_pathway=True):
        """Private function that returns the single parent produced species

        :param only_central: Focus on the central species only

        :type only_central: bool

        :return: List of node ids
        :rtype: list
        """
        return self.__isolatedSpecies(
            species=species,
            only_central=only_central,
            only_rp_pathway=only_rp_pathway,
            func_eq_0=self.G.successors,
            func_gt_0=self.G.predecessors
        )


    ## Recursive function that finds the order of the reactions in the graph
    #
    # NOTE: only works for linear pathways... need to find a better way ie. Tree's
    #
    def _recursiveReacSuccessors(self, node_name, reac_list, all_res, num_reactions):
        current_reac_list = [i for i in reac_list]
        self.logger.debug('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        succ_node_list = [i for i in self.G.successors(node_name)]
        flat_reac_list = [i for sublist in reac_list for i in sublist]
        self.logger.debug('flat_reac_list: '+str(flat_reac_list))
        self.logger.debug('current_reac_list: '+str(current_reac_list))
        if len(flat_reac_list)==num_reactions:
            self.logger.debug('Returning')
            #return current_reac_list
            all_res.append(current_reac_list)
            return all_res
        self.logger.debug('succ_node_list: '+str(succ_node_list))
        if not succ_node_list==[]:
            #can be multiple reactions at a given step
            multi_reac = []
            for n_n in succ_node_list:
                n = self.G.nodes.get(n_n)
                if n['type']=='reaction':
                    if not n_n in flat_reac_list:
                        multi_reac.append(n_n)
            #remove the ones that already exist
            self.logger.debug('multi_reac: '+str(multi_reac))
            multi_reac = [x for x in multi_reac if x not in flat_reac_list]
            self.logger.debug('multi_reac: '+str(multi_reac))
            if multi_reac:
                current_reac_list.append(multi_reac)
            self.logger.debug('current_reac_list: '+str(current_reac_list))
            #loop through all the possibilities
            for n_n in succ_node_list:
                n = self.G.nodes.get(n_n)
                if n['type']=='reaction':
                    if n_n in multi_reac:
                        self._recursiveReacSuccessors(n_n, current_reac_list, all_res, num_reactions)
                elif n['type']=='species':
                    if n['rp_trunk_species']==True:
                        self._recursiveReacSuccessors(n_n, current_reac_list, all_res, num_reactions)
        return all_res


    ##
    #
    # NOTE: only works for linear pathways... need to find a better way
    #
    def _recursiveReacPredecessors(self, node_name, reac_list, all_res, num_reactions):
        current_reac_list = [i for i in reac_list]
        self.logger.debug('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        pred_node_list = [i for i in self.G.predecessors(node_name)]
        flat_reac_list = [i for sublist in reac_list for i in sublist]
        self.logger.debug('flat_reac_list: '+str(flat_reac_list))
        self.logger.debug('current_reac_list: '+str(current_reac_list))
        if len(flat_reac_list)==num_reactions:
            self.logger.debug('Returning')
            #return current_reac_list
            all_res.append(current_reac_list)
            return all_res

            #can be multiple reactions at a given step
            multi_reac = []
            for n_n in pred_node_list:
                n = self.G.nodes.get(n_n)
                if n['type']=='reaction':
                    if not n_n in flat_reac_list:
                        multi_reac.append(n_n)
            #remove the ones that already exist
            self.logger.debug('multi_reac: '+str(multi_reac))
            multi_reac = [x for x in multi_reac if x not in flat_reac_list]
            self.logger.debug('multi_reac: '+str(multi_reac))
            if multi_reac:
                current_reac_list.append(multi_reac)
            self.logger.debug('current_reac_list: '+str(current_reac_list))
            #loop through all the possibilities
            for n_n in pred_node_list:
                n = self.G.nodes.get(n_n)
                if n['type']=='reaction':
                    if n_n in multi_reac:
                        self._recursiveReacPredecessors(n_n, current_reac_list, all_res, num_reactions)
                elif n['type']=='species':
                    if n['rp_trunk_species']==True:
                        self._recursiveReacPredecessors(n_n, current_reac_list, all_res, num_reactions)
        return all_res


    '''
    def _recursiveHierarchy(self, node_name, num_nodes, ranked_nodes):
        self.G.successors(node_name)
    '''

    ######################################################################################################
    ########################################## Public Function ###########################################
    ######################################################################################################


    def _recursiveReacPredecessors(self, node_name, reac_list):
        """Return the next linear predecessors

        Better than before, however bases itself on the fact that central species do not have multiple predesessors
        if that is the case then the algorithm will return badly ordered reactions

        :param node_name: The id of the starting node
        :param reac_list: The list of reactions that already have been run

        :type node_name: str
        :type reac_list: list

        :return: List of node ids
        :rtype: list
        """
        self.logger.debug('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        pred_node_list = [i for i in self.G.predecessors(node_name)]
        self.logger.debug(pred_node_list)
        if pred_node_list==[]:
            return reac_list
        for n_n in pred_node_list:
            n = self.G.nodes.get(n_n)
            if n['type']=='reaction':
                if n_n in reac_list:
                    return reac_list
                else:
                    reac_list.append(n_n)
                    self._recursiveReacPredecessors(n_n, reac_list)
            elif n['type']=='species':
                if n['rp_trunk_species']==True:
                    self._recursiveReacPredecessors(n_n, reac_list)
                else:
                    return reac_list
        return reac_list


    def orderedRetroReactions(self):
        """Public function to return the linear list of reactions

        :return: List of node ids
        :rtype: list
        """
        #Note: may be better to loop tho
        for prod_spe in self._onlyProducedSpecies():
            self.logger.debug('Testing '+str(prod_spe))
            ordered = self._recursiveReacPredecessors(prod_spe, [])
            self.logger.debug(ordered)
            if len(ordered)==self.num_reactions:
                return [i for i in reversed(ordered)]
        self.logger.error('Could not find the full ordered reactions')
        return []

    ############################# graph analysis ################################


    def exportJSON(self):
        return json_graph.node_link_data(self.G)


    ## Warning that this search algorithm only works for mono-component that are not networks (i.e where reactions follow each other)
    # DEPRECATED: this is linear
    # NOTE: only works for linear pathways... need to find a better way
    #
    def orderedRetroReactions(self):
        #Note: may be better to loop tho
        succ_res = []
        for cons_cent_spe in self._onlyConsumedSpecies():
            res = self._recursiveReacSuccessors(cons_cent_spe, [], [], self.num_reactions)
            if res:
                self.logger.debug(res)
                if len(res)==1:
                    succ_res = res[0]
                else:
                    self.logger.error('Multiple successors results: '+str(res))
            else:
                self.logger.warning('Successors no results')
        prod_res = []
        for prod_cent_spe in self._onlyProducedSpecies():
            res = self._recursiveReacPredecessors(prod_cent_spe, [], [], self.num_reactions)
            if res:
                self.logger.debug(res)
                if len(res)==1:
                    prod_res = [i for i in reversed(res[0])]
                else:
                    self.logger.error('Mutliple predecessors results: '+str(res))
            else:
                self.logger.warning('Predecessors no results')
        if succ_res and prod_res:
            if not succ_res==prod_res:
                self.logger.warning('Both produce results and are not the same')
                self.logger.warning('succ_res: '+str(succ_res))
                self.logger.warning('prod_res: '+str(prod_res))
            else:
                self.logger.debug('Found solution: '+str(succ_res))
                return succ_res
        return []
