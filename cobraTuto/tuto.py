#
#
# Maxime Lecomte
# 17/ 02 / 2020
# first script using cobra
#
#

from __future__ import print_function
import cobra
import cobra.test
import os
from os.path import join
from cobra.test import test_all
from cobra.test import create_test_model
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis



def test():
    """
    Test function.

    Returns
    -------
    type int
    0 if correctly installed.
    'exit code.OK'
    """
    print(test_all())

class TutoTest():
    """docstring forTutoTest."""

    def tuto(self):
        """
        tutorial cobra py

        """
        cobra_config = cobra.Configuration() # get configuration object
        print(cobra_config.lower_bound)
        print(cobra_config.upper_bound)
        print(cobra_config.bounds)

        cobra_config.bounds = -10, 20
        print(cobra.Reaction("R1"))
        print(cobra_config.bounds)
        print(cobra.Reaction("R2", lower_bound=None)) # set reversible reaction

        model = create_test_model("textbook")

        print(model.reactions.ACt2r)
        model.solver ## optlang.glpk_interface.Model object at 0x10110d3c8
        cobra_config.solver = "glpk_exact"
        new_model = create_test_model("textbook")
        new_model.solver ## optlang.glpk_exact_interface.Model object at 0x126bc75c0

    def addToModel(self):
        """
        this function allows to :
            1. add metabolites to a reaction
            2. add reaction to a model
        """

        # Best practise: SBML compliant IDs
        model = Model('example_model')

        reaction = Reaction('3OAS140')
        reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
        reaction.subsystem = 'Cell Envelope Biosynthesis'
        reaction.lower_bound = 0.  # This is the default
        reaction.upper_bound = 1000.  # This is the default

        ##### if i want a sepcific metabolit of model i can use : model.get_by_id

        ACP_c = Metabolite('ACP_c',formula='C11H21N2O7PRS',name='acyl-carrier-protein',compartment='c')
        omrsACP_c = Metabolite('3omrsACP_c',formula='C25H45N2O9PRS',name='3-Oxotetradecanoyl-acyl-carrier-protein', compartment='c')
        co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
        malACP_c = Metabolite('malACP_c',formula='C14H22N2O10PRS',name='Malonyl-acyl-carrier-protein',compartment='c')
        h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
        ddcaACP_c = Metabolite('ddcaACP_c',formula='C23H43N2O8PRS',name='Dodecanoyl-ACP-n-C120ACP',compartment='c')

        reaction.add_metabolites({
        malACP_c: -1.0,
        h_c: -1.0,
        ddcaACP_c: -1.0,
        co2_c: 1.0,
        ACP_c: 1.0,
        omrsACP_c: 1.0
        })

        print(reaction.name)
        reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
        print(reaction.genes)

        print('\n ############## without adding to a model #############\n')
        print('%i reactions initially' % len(model.reactions))
        print('%i metabolites initially' % len(model.metabolites))
        print('%i genes initially' % len(model.genes))
        ## no one in model so add it

        print("\n ###############  after adding to a model ###############\n")
        model.add_reactions([reaction])
        print('%i reaction' % len(model.reactions))
        print('%i metabolites' % len(model.metabolites))
        print('%i genes' % len(model.genes))
        ## is not empty

        print("\nReactions")
        print("---------")

        print("\nGenes and reaction linked to gene \n")

        for react in model.reactions:
            print(str(react.id) + ' : ' + str(react.reaction))
            # print("%s : %s" % (x.id, x.reaction)) ## working

            for gene in model.genes:
                print(str(gene.id) + ' : ' + str(react.id))

        print("\nMetabolites")
        print("---------")
        for met in model.metabolites:
            print(str(met.id) + ' : ' + str(met.formula))

        print('\n')
        model.objective = '3OAS140'
        print(model.objective.expression)
        print(model.objective.direction) ## maximize the flux in forward direction

    def readingWritting(self):
        data_dir = cobra.test.data_dir
        print("mini test files: ")
        print(", ".join(i for i in os.listdir(data_dir) if i.startswith("mini")))
        textbook_model = cobra.test.create_test_model("textbook")
        ecoli_model = cobra.test.create_test_model("ecoli")
        salmonella_model = cobra.test.create_test_model("salmonella")

        cobra.io.read_sbml_model(join(data_dir, "mini_fbc2.xml")) #mini_textbook
        cobra.io.write_sbml_model(textbook_model, "test_fbc2.xml") # generate a valid SBML file according id of the model


        # cobra.io.read_sbml_model(join(data_dir, "mini_cobra.xml"))
        # cobra.io.write_sbml_model(textbook_model, "test_cobra.xml",use_fbc_package=False) ## working with no flag use_fbc_package


        cobra.io.load_json_model(join(data_dir, "mini.json"))
        cobra.io.save_json_model(textbook_model, "test.json")


        cobra.io.load_yaml_model(join(data_dir, "mini.yml"))
        cobra.io.save_yaml_model(textbook_model, "test.yml")

        cobra.io.load_matlab_model(join(data_dir, "mini.mat"),variable_name="mini_textbook")
        cobra.io.save_matlab_model(textbook_model, "test.mat")

    def simulatingFBA(self):
        model = cobra.test.create_test_model("textbook")
        solution = model.optimize()
        # print(solution.status)
        # print(solution.fluxes)
        # print(solution.shadow_prices)
        # print(model.summary())
        # print(model.metabolites.nadh_c.summary())
        # print(model.metabolites.atp_c.summary())
        biomass_rxn = model.reactions.get_by_id('Biomass_Ecoli_core')
        # print(linear_reaction_coefficients(model)) ## 1
        # change the objective to ATPM
        model.objective = "ATPM"

        # The upper bound should be 1000, so that we get
        # the actual optimal value
        model.reactions.get_by_id("ATPM").upper_bound = 1000.
        # print(linear_reaction_coefficients(model))
        # print(model.summary()) ## objective id has changed it's ATPM

        print(flux_variability_analysis(model, model.reactions[:10]))







if __name__ == '__main__':
    # test() ## if you want to test if all is well installed
    # TutoTest().tuto()
    # TutoTest().addToModel()
    # TutoTest().readingWritting()
    TutoTest().simulatingFBA()
