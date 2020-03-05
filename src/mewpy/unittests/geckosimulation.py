"""
##################################################################

GECKO over REFRAMED tests

##################################################################
"""
from mewpy.model.gecko import GeckoModel, ModelList
from mewpy.simulation.reframed import GeckoSimulation
from mewpy.simulation import SimulationMethod
from collections import OrderedDict
import statistics
import pandas as pd
import os

def simulation_one():
    """
    Identify all pairs of (reaction,reaction_REV) associated to each protein draw reaction
    """
    model = GeckoModel('single-pool')
    simulation = GeckoSimulation(model)
    
    result= simulation.simulate(method= SimulationMethod.pFBA)
    wt_fluxes = result.fluxes

    rev_pairs = model.protein_rev_reactions
    with open("protein_reaction_under.csv",'w') as f:
        f.write("rxn; rxn_REV ; ; WT_rxn_flux; WT_rxn_REV_flux; ; O_rxn_flux;O_rxn_REV_flux\n")
        for protein_id in rev_pairs.keys():
            constraints = OrderedDict()
            rxn ='draw_prot_{}'.format(protein_id)
            constraints[rxn]= (0,0.5*wt_fluxes[rxn]) 
            result= simulation.simulate(constraints= constraints, method= SimulationMethod.pFBA)
            ssfluxes = result.fluxes
            if ssfluxes:
                f.write(protein_id+";;;;;;\n")
                for r_id, r_rev_id in rev_pairs[protein_id]:    
                    f.write("{}; {}; ; {}; {}; ;{}; {}\n".format(r_id,r_rev_id,wt_fluxes[r_id],wt_fluxes[r_rev_id],ssfluxes[r_id],ssfluxes[r_rev_id]))


def simulation_two():
    """
    Saves essential proteins to file
    """
    model = GeckoModel('single-pool')
    simulation = GeckoSimulation(model)
    essential_prot = simulation.essential_proteins
    with open('target-prot-single-pool.txt','w') as f:
        for p in model.proteins:
            if p not in essential_prot:
                f.write(p+"\n")




def simulation_three():
    constraints ={'draw_prot_Q04396': (2.3029771980628338e-07,1000),'draw_prot_P27515':(0,0.0),'draw_prot_P53687':(0.0,1000),'draw_prot_P81450': 0,
                  'draw_prot_P41939': (3.944331158229099e-06,1000),	'r_0659_REVNo1':0,'draw_prot_P19657':(0,0.0),'draw_prot_P37254':(0.0,1000),
                  'draw_prot_P33317': (0,0.0),'draw_prot_P54114':(0,0.0),'draw_prot_P25373':(0.0,1000),	'draw_prot_P08067':(0,2.1621430284827036e-06),
                  'draw_prot_P39002': (0.0,1000),'draw_prot_Q04792':(0.0,1000),	'draw_prot_Q12320':	(0.0,1000),	'draw_prot_P07347':(0.0,1000), 'draw_prot_P00958':(0,3.4862363761891556e-06),
                  'draw_prot_P15179': (0,0.0), 'draw_prot_P53332':(0.0,1000),'draw_prot_P23202':(0.0,1000),	'draw_prot_P21672':	(0,0.0),'draw_prot_P32191':(0,0.0),
                  'draw_prot_P38071': 0,'draw_prot_P06169':	(0.00019588347980199452,1000),'draw_prot_P22803':0,	'draw_prot_P38986':(0.0,1000),
                  'draw_prot_P32621': (0.0,1000),'draw_prot_P27472':(1.7311263946838522e-07,1000),'draw_prot_P00549': (0,1.624598159388334e-05),
                  'draw_prot_P32796':0}
    model = GeckoModel('single-pool')
    model.set_objective({'r_2111': 0.0 , 'r_4041':1.0})
    simulation = GeckoSimulation(model)
    result = simulation.simulate(method=SimulationMethod.pFBA)
    reference = result.fluxes
    
    result = simulation.simulate(constraints=constraints)
    
    from mewpy.optimization.evaluation import WYIELD, BPCY, TargetFlux
    evaluator_1 = WYIELD("r_2111", "r_2056", parsimonious = True)
    print(evaluator_1.get_fitness(result,None))
    evaluator_2 = BPCY("r_2111", "r_2056", "r_1714_REV", method=SimulationMethod.lMOMA ,reference=reference)
    print(evaluator_2.get_fitness(result,None))
    evaluator_3 = TargetFlux("r_2056")
    print(evaluator_3.get_fitness(result,None))


def test_basic_gecko_adjustment():
    in_model = {'P00549': 0.1, 'P31373': 0.1, 'P31382': 0.1, 'P39708': 0.1, 'P39714': 0.1, 'P39726': 0.1, 'Q01574': 0.1}
    not_in_model = {'P10591': 0.1, 'P31383': 0.1, 'P32471': 0.1}
    measurements = pd.concat([pd.Series(in_model), pd.Series(not_in_model)])
    model = GeckoModel('multi-pool')
    model.limit_proteins(fractions=pd.Series(measurements))
    simul = GeckoSimulation(model)
    sol = simul.simulate()
    assert sol.objective_value > 0.05
    assert len(model.proteins) - len(model.pool_proteins) - len(in_model) == 0
    assert all(model.reactions[rxn].ub > 0 for rxn in model.individual_protein_exchanges)




def test_gecko_adjustment_sanchez_etal():
    mmol_gdw = pd.read_csv(os.path.join(os.path.dirname(__file__), '../model/data/sanchez-mmol_gdw.csv'))
    PROTEIN_PROPERTIES = ModelList().protein_properties()
    ggdw = pd.Series(PROTEIN_PROPERTIES.loc[mmol_gdw.index, 'mw'] / 1000.) * pd.Series(mmol_gdw)
    model = GeckoModel('multi-pool')
    simulation = GeckoSimulation(model)
    result = simulation.simulate()
    growth_rate_unlimited_protein = result.objective_value
    model.limit_proteins(ggdw=pd.Series(ggdw))
    result = simulation.simulate()
    growth_rate_limited_protein = result.objective_value
    # should be smaller, but how much..
    assert growth_rate_limited_protein < 0.8 * growth_rate_unlimited_protein
    measured_in_model = set(mmol_gdw.index).intersection(model.proteins)
    assert sum(model.concentrations[p] - ggdw[p] for p in measured_in_model) < 1e-10
    #assert sum(abs(rxn.upper_bound - mmol_gdw[rxn.metadata['uniprot']])
    #           for rxn in model.individual_protein_exchanges) < 1e-6
    #assert sum(rxn.metabolites[model.common_protein_pool] +
    #           PROTEIN_PROPERTIES.loc[rxn.annotation['uniprot'], 'mw'] / 1000.
    #           for rxn in model.pool_protein_exchanges) < 1e-6
    #assert model.p_measured > 0.25                                  # With yeast 8.1.3 -> p_measured = 0.296
    #assert model.f_mass_fraction_measured_matched_to_total > 0.25   # With yeast 8.1.3 -> f = 0.304
    #assert model.protein_pool_exchange.upper_bound > 0.015          # With yeast 8.1.3 -> pool_exchange = 0.0212



def simulation_four():

    constraints = {'draw_prot_P17505':0,
                   'draw_prot_P54885':0,
                   'draw_prot_P50107':0,
                   'draw_prot_P37303':0,		
                   'draw_prot_P38113':0,
                   'draw_prot_Q06408':0,
                   'draw_prot_P32340':0,
                   'draw_prot_P00817':0,
                   'draw_prot_P09440':0,
                   'draw_prot_P42951':0,
                   'draw_prot_P36013':0,
                   'draw_prot_P32473':0,
                   'draw_prot_P27680':0,		
                   'draw_prot_P41939':0,
                   'draw_prot_P17695':0,
                   'draw_prot_P32383':0,		
                   'draw_prot_P38840':0,
                   'draw_prot_P38715':0,
                   'draw_prot_P23542':0,		
                   'draw_prot_Q01574':0,		
                   'draw_prot_P06208':0,	
                   'draw_prot_P00330':0,	
                   'draw_prot_P33330':0,		
                   'draw_prot_P32419':0,		
                   'draw_prot_P47143':0,		
                   'draw_prot_P32179':0,	
                   'draw_prot_P06169':0}
    
    model = GeckoModel('single-pool')
    model.set_objective({'r_2111': 0 , 'r_4041':1.0})
    simulation = GeckoSimulation(model)
    result = simulation.simulate(method = SimulationMethod.pFBA)
    reference = result.fluxes
    print("biomass {}   tyrosine: {}".format(reference['r_4041'],reference['r_1913']))
    result = simulation.simulate(method = SimulationMethod.pFBA,constraints=constraints)
    reference = result.fluxes
    print("biomass {}   tyrosine: {}".format(reference['r_4041'],reference['r_1913']))

    reactions = [] 
    for prot_exchange in constraints.keys():
        prot = prot_exchange[len('draw_prot_'):]
        reactions.append(simulation.protein_reactions(prot))
    reactions.sort()    
    print(reactions)
    print(len(reactions))


if __name__ == "__main__":
    #simulation_one()
    #simulation_two()
    #simulation_three()
    #test_basic_gecko_adjustment()
    #test_gecko_adjustment_sanchez_etal()
    simulation_four()
    
    
