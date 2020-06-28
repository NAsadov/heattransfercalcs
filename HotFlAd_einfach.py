from math import exp, pi, log
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
import requests
from os.path import join

'''
Importing measurement data
'''
data_folder = '/Users/nasir/tubCloud/Shared/HotFlAd/03_Messdaten/Monitoring/data/8h'
temp_path = join(data_folder, 'arrTemps.txt')
MID_path = join(data_folder, 'arrMID.txt')
header_temp = ['bt'+str(i) for i in range(1, 20)]
header_MID = ['mid'+str(i) for i in range(1, 7)]
temp_data = pd.read_csv(temp_path, sep='\t', header=0, names=header_temp)
MID_data = pd.read_csv(MID_path, sep='\t', header=0, names=header_MID)
print(MID_data.head())


"""
Initialization of the temperature sensors

#t1 = 333.15  # [K] Temperature at the server HX inlet (pump box circuit)
#t2 = 332.15  # [K] Temperature at the server HX outlet (pump box circuit)

bt1 = 323.15  # [K] T at HT outlet (AdKM)
dt_1_5_meas = -0.5
bt5 = 322.65  # bt1 + dt_1_5_meas   # [K] T at HT outlet junction (real and sim server)
dt_5_17_meas = -0.5
bt17 = bt5 + dt_5_17_meas   # [K] Temperature at the buffer storage inlet
dt_17_18_meas = 12
bt18 = bt17 + dt_17_18_meas   # [K] Temperature at the buffer storage outlet
dt_5_t3_meas = -0.5
t3 = bt5 + dt_5_t3_meas   # [K] Temperature at the server HX inlet (HT circuit)
dt_t3_t4_meas = 12
t4 = t3 + dt_t3_t4_meas   # [K] Temperature at the server HX outlet (HT circuit)
dt_t4_6_meas = -0.5
bt6 = t4 + dt_t4_6_meas  # [K] T at HT inlet junction (real and sim server)
dt_6_2_meas = -0.5
bt2 = bt6 + dt_6_2_meas

# bt2 = 328.15  # [K] T at HT inlet (AdKM)

bt7 = 349.15  # [K] T at HX server-sim inlet (central heating circuit)
bt8 = 339.15  # [K] T at HX server-sim outlet (central heating circuit)
bt4 = 330.15  # [K] T at HX server-sim outlet (HT circuit)
bt3 = 322.15  # [K] T at HX server-sim inlet (HT circuit)

bt15 = 285.15  # [K] T re-cooling water in
bt16 = 292.15  # [K] T re-cooling water out
bt13 = 298.65  # [K] T re-cooled water in = T emulator mt out
bt14 = 293.15  # [K] T re-cooled water out
bt12 = 300.15  # [K] T emulator mt in

bt0 = 291.15  # [K] T at LT inlet (AdKM)
bt9 = 288.15  # [K] T at LT outlet (AdKM)
bt11 = 288.65  # [K] T at LT inlet to KAE
bt10 = 290.65  # [K] T at LT outlet from KAE
"""

"""
Initiation of the volume flow rate sensors

V_MID6 = 0.000167  # [m^3/s] Volume flow rate in the server circuit (10 l/min)
V_MID5 = 0.000167  # [m^3/s] Volume flow rate in the server secondary circuit (10 l/min)
V_MID4 = 0.000833  # [m^3/s] Volume flow rate in the server simulation circuit (50 l/min)
V_MID3 = 0.000833  # [m^3/s] Volume flow rate in the server simulation circuit (50 l/min) (= V_BF1 ?)
# V_BF1 = 0.000833  # [m^3/s] Volume flow rate at the adsorption chiller outlet (50 l/min)  # V_BF sensors not readable?
# V_BF4 = 0.000833  # [m^3/s] Volume flow rate of the re-cooling water (50 l/min)
V_MID2 = 0.001  # [m^3/s] Volume flow rate in the re-cooling circuit (60 l/min)
V_MID1 = 0.000833  # [m^3/s] Volume flow rate in the cold water circuit (50 l/min)
"""


dict_pipe_sizes = {25: {'ri': 0.0136, 'ro': 0.01685},  # DIN EN 10255
                   32: {'ri': 0.01795, 'ro': 0.0212},
                   40: {'ri': 0.0209, 'ro': 0.02415}}


'''class ThermodynamicProperties:
    def __init__(self, t, fluid='IF97::Water', p=101325):
        """
        This class calculates the cp and rho of the fluid from its temperature.

        :param fluid: which fluid
        :param p: pressure of the fluid [Pa]
        """
        self.t = t
        def get_cp(self, t):
            return PropsSI('C', 'T', self.t, 'P', p, fluid)

        def get_rho(self, t):
            return PropsSI('D', 'T', self.t, 'P', p, fluid)
'''

'''
class Pipe(ThermodynamicProperties):
    def __init__(self, t_in, dn, r3, material, t_amb=295.15):
        """
        This class calculates the output temperature of the pipe.

        :param t_in: inlet temperature [K]
        :param t_amb: ambient temperature [K]
        :param r1: inner radius of the pipe [m]
        :param r2: outer radius of the pipe [m]
        :param r3: outer radius of the insulated pipe [m]
        :param h1: heat transfer coefficient between the fluid and the inner surface of the pipe [W/m2K]
        :param h2: heat transfer coefficient between the outer surface of the insulated pipe and ambient [W/m2K]
        :param k1: thermal conductivity of the pipe material [W/mK]
        :param k2: thermal conductivity of the insulation material [W/mK]
        """
        super().__init__(fluid, p)
        self.t_in = t_in
        self.t_amb = t_amb
        self.dn = dn
        self.r1 = dict_pipe_sizes[dn]['ri']
        self.r2 = dict_pipe_sizes[dn]['ro']
        self.r3 = self.r2 + 0.009
        self.k1 = k1
        self.k2 = k2
        
        self.A1 = 2 * pi * self.r1 * self.L
        self.A3 = 2 * pi * self.r3 * self.L
        
        self.Ri = 1 / (self.h1 * self.A1)
        self.R1 = log(self.r2 / self.r1) / (2 * pi * self.k1 * self.L)
        self.R2 = log(self.r3 / self.r2) / (2 * pi * self.k2 * self.L)
        self.Ro = 1 / (self.h2 * self.A3)
        self.R_tot = np.sum(self.Ri, self.R1, self.R2, self.Ro)

        self.t_out = self.t_in
        for dL in np.linspace(0, self.L, 10):
            self.Q_loss = (self.t_out - self.t_amb) / self.R_tot
            self.t_out = self.t_in - self.Q_loss * dL / (get_cp(self.t_out) * get_rho(self.t_out) * self.A1 * dL)
'''


class Node:
    def __init__(self, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in, p=101325, fluid='IF97::Water'):
        """
        This class calculates the heat transfer via measured values.

        :param tc_in_meas: Measured temperature at the colder circuit inlet [K]
        :param tc_out_meas: Measured outlet temperature at colder circuit
        :param dvc_in:
        :param th_in_meas:
        :param th_out_meas:
        :param dvh_in:
        :param p:
        :param fluid:
        """
        self.tc_in_meas = tc_in_meas
        self.tc_out_meas = tc_out_meas
        self.dtc_meas = self.tc_out_meas - self.tc_in_meas
        self.dVc_in = dvc_in
        self.cpc_in = PropsSI('C', 'T', self.tc_in_meas, 'P', p, fluid)
        self.rhoc_in = PropsSI('D', 'T', self.tc_in_meas, 'P', p, fluid)
        self.dmc_in = self.rhoc_in * self.dVc_in
        self.Qc_meas = round(self.dmc_in * self.cpc_in * self.dtc_meas)

        self.th_in_meas = th_in_meas
        self.th_out_meas = th_out_meas
        self.dth_meas = self.th_in_meas - self.th_out_meas
        self.dVh_in = dvh_in
        self.cph_in = PropsSI('C', 'T', self.th_in_meas, 'P', p, fluid)
        self.rhoh_in = PropsSI('D', 'T', self.th_in_meas, 'P', p, fluid)
        self.dmh_in = self.rhoh_in * self.dVh_in
        self.Qh_meas = round(self.dmh_in * self.cph_in * self.dth_meas)


class EpsNtu(Node):
    def __init__(self, ua, tc_in_meas, tc_out_meas, dvc_in, cpc_in, th_in_meas, th_out_meas, dvh_in, cph_in):
        """
        This class calculates heat transfer and the output temperatures via epsilon-NTU method.

        :param ua: UA value of the heat exchanger [W/K]
        :param tc_in_meas: Measured temperature at the colder circuit inlet [T]
        :param tc_out_meas:
        :param dvc_in:
        :param th_in_meas:
        :param th_out_meas:
        :param dvh_in:
        """
        super().__init__(tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in, cpc_in, cph_in)
        self.UA = ua
        self.dt_in = self.th_in_meas - self.tc_in_meas
        self.Cc = self.cpc_in * self.dmc_in
        self.Ch = self.cph_in * self.dmh_in
        self.Cmin = min(self.Cc, self.Ch)
        self.Cmax = max(self.Cc, self.Ch)
        self.Cr = self.Cmin / self.Cmax
        self.Qmax = self.Cmin * self.dt_in
        self.NTU = self.UA / self.Cmin
        self.eps = (1 - exp(-self.NTU * (1 - self.Cr))) / (1 - self.Cr * exp(-self.NTU * (1 - self.Cr)))
        self.Q_sim = self.Qmax * self.eps
        self.tc_out_sim = self.Q_sim / (self.dmc_in * self.cpc_in) + self.tc_in_meas
        self.th_out_sim = -self.Q_sim / (self.dmc_in * self.cph_in) + self.th_in_meas


class TableResults(EpsNtu):
    def __init__(self, ua, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in):
        """
        This class creates a dictionary and subsequently a pandas data frame with results.

        :param ua:
        :param tc_in_meas:
        :param tc_out_meas:
        :param dvc_in:
        :param th_in_meas:
        :param th_out_meas:
        :param dvh_in:
        """
        super().__init__(ua, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in)
        self.results_dict = {'T inlet': [self.tc_in_meas, self.th_in_meas],
                             'T outlet': [self.tc_out_meas, self.th_out_meas],
                             'Delta T': [self.dtc_meas, self.dth_meas],
                             'Vol. fl. r. inlet': [self.dVc_in, self.dVh_in],
                             'Cp inlet': [self.cpc_in, self.cph_in],
                             'Rho inlet': [self.rhoc_in, self.rhoh_in],
                             'Mass fl. r.': [self.dmc_in, self.dmh_in],
                             'Simulation UA value': [self.UA, self.UA],
                             'Measured heat transfer': [self.Qc_meas, self.Qh_meas],
                             'Simulated heat transfer': [self.Q_sim, self.Q_sim]}
        self.columns = list(self.results_dict.keys())
        self.index = ['Cold side', 'Hot side']
        self.results_df = pd.DataFrame(self.results_dict, columns=self.columns, index=self.index)


"""
Sandbox
"""
emul_df = TableResults(650, temp_data.bt11, temp_data.bt10, MID_data.mid1, temp_data.bt12, temp_data.bt13, MID_data.mid2)
print(emul_df.results_df)
#serv_df = TableResults(3250, t3, t4, mid6, t1, t2, mid5)
#print(serv_df.results_df)
#recooler_df = TableResults(3500, bt15, bt16, V_BF4, bt13, bt14, mid2)
#print(recooler_df.results_df)
