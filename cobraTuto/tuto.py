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
import pandas as pd
import pylab
import matplotlib
import matplotlib.pyplot as plt
import numpy as np



from os.path import join
from cobra.test import test_all
from cobra.test import create_test_model
from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis
from time import time
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
# from matplotlib import plot_helper




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

        print(textbook_model.summary())
        print(ecoli_model.summary())
        # print(salmonella_model.summary()) ### doresn't work

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
        print(flux_variability_analysis(model, model.reactions[:10], fraction_of_optimum=0.9))

        loop_reactions = [model.reactions.FRD7, model.reactions.SUCDi]
        flux_variability_analysis(model, reaction_list=loop_reactions, loopless=False)   ## see chapter 5 on optimizing methods in metabolic methods textbook

        flux_variability_analysis(model, reaction_list=loop_reactions, loopless=True)  ## FRD7 can be remove

        model.optimize()
        print(model.summary(fva=0.95))  ## metabolite mass balance can be also check with FVA
        # print(model.summary())

        print(model.metabolites.pyr_c.summary(fva=0.95))  # see for pyr_c, where he is using in product ( producing), and use as a substrate ( consumming) in some reactions ( rxn id)

        ##### parsimonious FBA #####

        model.objective = 'Biomass_Ecoli_core'
        fba_solution = model.optimize()
        pfba_solution = cobra.flux_analysis.pfba(model)
        # print(pfba_solution.objective_value)
        # print(pfba_solution.status)
        # print(pfba_solution.fluxes)
        # print(pfba_solution.shadow_prices)

        ##### geometric FBA finds a unique optimal flux distribution which is central to the range of possible fluxes

        geometric_fba_sol = cobra.flux_analysis.geometric_fba(model)
        # geometric_fba_sol

    def simulatingDeletions(self):
        cobra_model = cobra.test.create_test_model("textbook")
        ecoli_model = cobra.test.create_test_model("ecoli")

        #### knocking out single genes and reactions #####
        print('complete model: ', cobra_model.optimize())
        with cobra_model:
            cobra_model.reactions.PFK.knock_out()
            print('pfk knocked out: ',cobra_model.optimize())
            print(cobra_model.optimize().status)  ## optimal

        print('complete model: ', cobra_model.optimize())
        with cobra_model:
            cobra_model.reactions.PFK.knock_out()
            print('pfk knocked out: ', cobra_model.optimize())
            print(cobra_model.optimize().status)

        print('complete model: ', cobra_model.optimize())
        with cobra_model:
            cobra_model.genes.b1723.knock_out()
            print('pfkA knocked out: ', cobra_model.optimize())
            print(cobra_model.optimize().status)
            cobra_model.genes.b3916.knock_out()
            print('pfkB knocked out: ', cobra_model.optimize())
            print(cobra_model.optimize().status)

        ##### single deletions ####
        deletion_result = single_gene_deletion(cobra_model, cobra_model.genes[:20]) ## also be done for all or subset reaction, all/subset genes
        print(deletion_result)

        ##### Double Deletions ####

        double_deletion_result  = double_gene_deletion(
        cobra_model, cobra_model.genes[:20]).round(4)
        print(double_deletion_result)

        start = time()  # start timer()
        double_gene_deletion(
            ecoli_model, ecoli_model.genes[:25], processes=2)
        t1 = time() - start
        print("Double gene deletions for 200 genes completed in "
              "%.2f sec with 2 cores" % t1)

        start = time()  # start timer()
        double_gene_deletion(
            ecoli_model, ecoli_model.genes[:25], processes=1)
        t2 = time() - start
        print("Double gene deletions for 200 genes completed in "
              "%.2f sec with 1 core" % t2)

        print("Speedup of %.2fx" % (t2 / t1))

        ## double reaction can also be run for reaction it's the same thing

    def productionEnvelops(self):  ## allow to show distinct pahses of optimal growth with different use of two different substrates
        model = cobra.test.create_test_model("textbook")
        prod_env = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
        prod_env.head()
        prod_env = production_envelope(
            model, ["EX_o2_e"], objective="EX_ac_e", carbon_sources="EX_glc__D_e")
        prod_env.head()
        prod_env.plot(kind='line',x='EX_o2_e',y='carbon_yield_maximum')
        plt.show()

    def fluxSampling(self):
        print("doesn't work cause of sample package")

    def looplessSolution(self):
        ''' identification of a thermodybamically consistent flux state whitout loops'''


        coli = cobra.test.create_test_model('ecoli') ## do not take salmonella model --> doesn't work

        nominal = coli.optimize()
        loopless = loopless_solution(coli) # identify the closest flux distribution  to the original one
        # print(loopless.status) # 0 and optimal

        df=pd.DataFrame(dict( loopless=loopless.fluxes,nominal=nominal.fluxes))
        plt.plot("loopless", "nominal", data=df, linestyle='none', marker='o')
        plt.show()




        #### loopless model ####
         #use it if :
            # 1- combine a non linear objective with the loopless condition
            # 2- force the model to be infeasible in the presence of loops indépendent of the set reaction bounds

        #documentation_builder.plot_helper.plot_loop() ### probleme

        model = Model()
        model.add_metabolites([Metabolite(i) for i in "ABC"])
        model.add_reactions([Reaction(i) for i in ["EX_A", "DM_C", "v1", "v2", "v3"]])

        model.reactions.EX_A.add_metabolites({"A": 1})
        model.reactions.DM_C.add_metabolites({"C": -1})

        model.reactions.v1.add_metabolites({"A": -1, "B": 1})
        model.reactions.v2.add_metabolites({"B": -1, "C": 1})
        model.reactions.v3.add_metabolites({"C": -1, "A": 1})

        model.objective = 'DM_C'
        with model:
            add_loopless(model)
            solution = model.optimize()
        print("loopless solution: status = " + solution.status)
        print("loopless solution flux: v3 = %.1f" %solution.fluxes["v3"])

        solution = pfba(model)
        print("parsimonious solution: status = " + solution.status)
        print("loopless solution flux: v3 = %.1f" % solution.fluxes["v3"])

        model.reactions.v3.lower_bound = 1
        with model:
            add_loopless(model)
            try:
                solution = model.optimize()
            except:
                print('model is infeasible')

        solution = pfba(model)
        print("parsimonious solution: status = " + solution.status)
        print("loopless solution flux: v3 = %.1f" % solution.fluxes["v3"])

    def consistencyTesting(self):
        ''' can detect all the blocked reactions and also give us consistent networks '''
        test_model = cobra.Model("test_model")
        v1 = cobra.Reaction("v1")
        v2 = cobra.Reaction("v2")
        v3 = cobra.Reaction("v3")
        v4 = cobra.Reaction("v4")
        v5 = cobra.Reaction("v5")
        v6 = cobra.Reaction("v6")
        test_model.add_reactions([v1, v2, v3, v4, v5, v6])

        v1.reaction = "-> 2 A"
        v2.reaction = "A <-> B"
        v3.reaction = "A -> D"
        v4.reaction = "A -> C"
        v5.reaction = "C -> D"
        v6.reaction = "D ->"

        v1.bounds = (0.0, 3.0)
        v2.bounds = (-3.0, 3.0)
        v3.bounds = (0.0, 3.0)
        v4.bounds = (0.0, 3.0)
        v5.bounds = (0.0, 3.0)
        v6.bounds = (0.0, 3.0)

        test_model.objective = v6
        print(cobra.flux_analysis.find_blocked_reactions(test_model))


        ### use FASTCC ###
            ## you can expect to efficiently obtain an accurate consistent network

        consistent_model = cobra.flux_analysis.fastcc(test_model)
        print(consistent_model.reactions)

        ### same result

    def gapFilling(self):
        pass








if __name__ == '__main__':
    # test() ## if you want to test if all is well installed
    test = TutoTest()
    # test.tuto()
    # test.addToModel()
    # test.readingWritting()
    # test.simulatingFBA()
    # test.simulatingDeletions()
    # test.productionEnvelops()
    # test.fluxSampling()
    # test.looplessSolution()
    # test.consistencyTesting()
    test.gapFilling()
