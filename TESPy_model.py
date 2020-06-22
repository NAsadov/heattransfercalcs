from tespy.networks import network
# create a network object with water as fluid
fluid_list = ['water']
my_plant = network(fluids=fluid_list)

my_plant.set_attr(T_unit='C', p_unit='bar', h_unit='kJ / kg')



from tespy.components import (
    sink, source, pipe, heat_exchanger)

soserv = source('Serverwärme Vorlauf')
siserv = sink('Serverwärme Rücklauf')

sork = source('Rückkühlung Vorlauf')
sirk = sink('Rückkühlung Rücklauf')

wtserv = heat_exchanger('Wärmetauscher')
wtserv.set_attr(pr1=0.98, pr2=0.98, ttd_u=5)
akm = adsorption_chiller('Adsorptionskältemaschine')
rke = heat_exchanger('Rückkühlemulator')
kae = heat_exchanger('Kälteabnahmeemulator')

lHTein = pipe('Leitung HT Eingang')
lHTaus = pipe('Leitung HT Ausgang')
lMTein = pipe('Leitung MT Eingang')
lMTaus1 = pipe('Leitung MT Ausgang nach KAE')
lMTaus2 = pipe('Leitung MT KAE nach RKE')
lNTein = pipe('Leitung NT Eingang')
lNTaus = pipe('Leitung NT Ausgang')

lHTein.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.025, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)

lHTaus.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.025, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)

lMTein.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.04, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)

lMTaus1.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.025, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)

lMTaus2.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.025, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)

lNTein.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.025, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)

lNTaus.set_attr(ks=0.0000015, # pipe's roughness in meters
L=2, # length in m
D=0.025, # diameter in m
kA=10, # area independent heat transfer coefficient kA in W/K
Tamb=20) # ambient temperature of the pipe environment (ground temperature)


from tespy.connections import connection

soserv_wtserv = connection(soserv, 'out1', wtserv, 'in1')
wtserv_siserv = connection(wtserv, 'out1', siserv, 'in1')

wtserv_lHTein = connection(wtserv, 'out2', lHTein, 'in1')
lHTein_akm = connection(lHTein, 'out1', akm, 'in1')
akm_lHTaus = connection(akm, 'out1', lHTaus, 'in1')
lHTaus_wtserv = connection(lHTaus, 'out1', wtserv, 'in2')

sork_rke = connection(sork, 'out1', rke, 'in1')
rke_sirk = connection(rke, 'out1', sirk, 'in1')

rke_lMTein = connection(rke, 'out2', lMTein, 'in1')
lMTein_akm = connection(lMTein, 'out1', akm, 'in2')
akm_lMTaus1 = connection(akm, 'out1', lMTaus1, 'in1')
lMTaus1_kae = connection(lMTaus1, 'out1', kae, 'in1')
kae_lMTaus2 = connection(kae, 'out1', lMTaus2, 'in1')
lMTaus2_rke = connection(lMTaus2, 'out1', rke, 'in2')

kae_lNTein = connection(kae, 'out2', lNTein, 'in1')
lNTein_akm = connection(lNTein, 'out1', akm, 'in3')
akm_lNTaus = connection(akm, 'out3', lNTaus, 'in1')
lNTaus_kae = connection(lNTaus, 'out1', kae, 'in2')


my_plant.add_conns(soserv_wtserv, wtserv_siserv, wtserv_lHTein, lHTein_akm, akm_lHTaus, lHTaus_wtserv,
                   sork_rke, rke_sirk, rke_lMTein, lMTein_akm, akm_lMTaus1, lMTaus1_kae, kae_lMTaus2, lMTaus2_rke,
                   kae_lNTein, lMTein_akm, akm_lNTaus, lNTaus_kae)


my_plant.set_attr(iterinfo=True)
my_plant.solve(mode='design')
my_plant.print_results()
