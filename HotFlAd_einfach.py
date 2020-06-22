from math import exp
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI

"""
Initiation of the temperature sensors
"""


t1 = 333.15  # [K] Temperature at the server HX inlet (pump box circuit)
t2 = 332.15  # [K] Temperature at the server HX outlet (pump box circuit)


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
Initiation of the volume flow rate sensors
"""


V_MID6 = 0.000167  # [m^3/s] Volume flow rate in the server circuit (10 l/min)
V_MID5 = 0.000167  # [m^3/s] Volume flow rate in the server secondary circuit (10 l/min)
V_MID4 = 0.000833  # [m^3/s] Volume flow rate in the server simulation circuit (50 l/min)
V_MID3 = 0.000833  # [m^3/s] Volume flow rate in the server simulation circuit (50 l/min) (= V_BF1 ?)
# V_BF1 = 0.000833  # [m^3/s] Volume flow rate at the adsorption chiller outlet (50 l/min)  # V_BF sensors not readable?
# V_BF4 = 0.000833  # [m^3/s] Volume flow rate of the re-cooling water (50 l/min)
V_MID2 = 0.001  # [m^3/s] Volume flow rate in the re-cooling circuit (60 l/min)
V_MID1 = 0.000833  # [m^3/s] Volume flow rate in the cold water circuit (50 l/min)


"""
Calculation of intermediate values
"""


class Pipe:
    def __init__(self, t_in, t_out, dt):
        self.t_in = t_in
        self.t_out = t_out
        self.dt = dt
        if self.dt is None:
            self.dt = self.t_in - self.t_out
        elif self.t_out is None:
            self.t_out = self.t_in - self.dt
        else:
            print('Input "None" in case a value is unknown')


class Node:
    def __init__(self, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in, p=101325, fluid='IF97::Water'):
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


"""
epsilon-NTU
"""


class EpsNtu(Node):
    def __init__(self, ua, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in):
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
pipe_1_5 = Pipe(bt1, None, 0.5)
print(pipe_1_5.t_out)

emul_df = TableResults(650, bt11, bt10, V_MID1, bt12, bt13, V_BF2)
print(emul_df.results_df)
serv_df = TableResults(3250, t3, t4, V_MID6, t1, t2, V_MID5)
print(serv_df.results_df)
recooler_df = TableResults(3500, bt15, bt16, V_BF4, bt13, bt14, V_BF2)
print(recooler_df.results_df)
