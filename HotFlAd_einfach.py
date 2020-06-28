from math import exp, pi, log
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI
import requests
import os
import matplotlib.pyplot as plt
import seaborn as sns
'''
Importing measurement data
'''
data_folder = '/Users/nasir/tubCloud/Shared/HotFlAd/03_Messdaten/Monitoring/data/8h'
temp_path = os.path.join(data_folder, 'arrTemps.txt')
mid_path = os.path.join(data_folder, 'arrMID.txt')
time_path = os.path.join(data_folder, 'arrTime.txt')
q_path = os.path.join(data_folder, 'arrQ.txt')

header_temp = ['bt'+str(i) for i in range(1, 20)]
header_mid = ['mid'+str(i) for i in range(1, 7)]
header_q = ['q'+str(i) for i in range(1, 7)]

temp_data = pd.read_csv(temp_path, sep='\t', header=0, names=header_temp)
mid_data = pd.read_csv(mid_path, sep='\t', header=0, names=header_mid)
time_data = pd.read_csv(time_path)
time_data['Europe/Berlin'] = pd.to_datetime(time_data['# Unix Time'], unit='s').astype('datetime64[ns, Europe/Berlin]')
q_data = pd.read_csv(q_path, sep='\t', header=0, names=header_q)
# print(temp_data.isnull().values.any())
# print(MID_data.isnull().values.any())
print(q_data.head())

'''
Plotting the input data
'''
f, ax = plt.subplots(figsize=(5, 6))
sns.set_style("whitegrid")
g = sns.lmplot(x=time_data['Europe/Berlin'], y=temp_data.bt1, data=temp_data, aspect=2)
g = (g.set_axis_labels("Time [s]", "Temperature [K]").set(xlim=(0, len(time_data['Europe/Berlin'])), ylim=(0, 100)))
plt.title("Temperature profiles")
plt.show(g)


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
        def c2k(t):
            return t + 273.15

        def lpm2cmps(v):
            return v * 1.6667e-5

        self.tc_in_meas = c2k(tc_in_meas)
        self.tc_out_meas = c2k(tc_out_meas)
        self.dtc_meas = self.tc_out_meas - self.tc_in_meas
        self.dVc_in = lpm2cmps(dvc_in)
        self.cpc_in = PropsSI('C', 'T', self.tc_in_meas, 'P', p, fluid)
        self.rhoc_in = PropsSI('D', 'T', self.tc_in_meas, 'P', p, fluid)
        self.dmc_in = self.rhoc_in * self.dVc_in
        self.Qc_meas = round(self.dmc_in * self.cpc_in * self.dtc_meas)

        self.th_in_meas = c2k(th_in_meas)
        self.th_out_meas = c2k(th_out_meas)
        self.dth_meas = self.th_in_meas - self.th_out_meas
        self.dVh_in = lpm2cmps(dvh_in)
        self.cph_in = PropsSI('C', 'T', self.th_in_meas, 'P', p, fluid)
        self.rhoh_in = PropsSI('D', 'T', self.th_in_meas, 'P', p, fluid)
        self.dmh_in = self.rhoh_in * self.dVh_in
        self.Qh_meas = round(self.dmh_in * self.cph_in * self.dth_meas)


class EpsNtu(Node):
    def __init__(self, ua, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in):
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
        super().__init__(tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in)
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
        self.results_dict = {'T inlet': [tc_in_meas, th_in_meas],
                             'T outlet': [tc_out_meas, th_out_meas],
                             'Delta T': [self.dtc_meas, self.dth_meas],
                             'Vol. fl. r. inlet': [dvc_in, dvh_in],
                             'Cp inlet': [self.cpc_in, self.cph_in],
                             'Rho inlet': [self.rhoc_in, self.rhoh_in],
                             'Mass fl. r.': [self.dmc_in, self.dmh_in],
                             'Simulation UA value': [self.UA, self.UA],
                             'Measured heat transfer': [self.Qc_meas, self.Qh_meas],
                             'Imported heat transfer': [q_data.q1[sr], q_data.q2[sr]],
                             'Simulated heat transfer': [self.Q_sim, self.Q_sim]}
        self.columns = list(self.results_dict.keys())
        self.index = ['Cold side', 'Hot side']
        self.results_df = pd.DataFrame(self.results_dict, columns=self.columns, index=self.index)


"""
Sandbox
"""
sr = 2000  # Sample row number
emul_df = TableResults(450, temp_data.bt10[sr], temp_data.bt11[sr], mid_data.mid1[sr], temp_data.bt12[sr],
                       temp_data.bt13[sr], mid_data.mid2[sr])
print(emul_df.results_df)
#serv_df = TableResults(3250, t3, t4, mid6, t1, t2, mid5)
#print(serv_df.results_df)
#recooler_df = TableResults(3500, bt15, bt16, V_BF4, bt13, bt14, mid2)
#print(recooler_df.results_df)
