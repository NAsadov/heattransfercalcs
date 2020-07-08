from math import exp, pi, log
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns
from thermodyn.functions import get_rmse, get_air_property_at1bar, get_k, c2k
from thermodyn.classes.epsntu import EpsNtu
from thermodyn.classes.pipe_heat_loss import Pipe

'''
Importing measurement data
'''
monitoring_folder = '/Users/nasir/tubCloud/Shared/HotFlAd/03_Messdaten/Monitoring'
simulation_folder = os.path.join(monitoring_folder, 'Simulation')
temporary_path = os.path.join(monitoring_folder, 'ReglerOptimierungKK', '2nd_testing', 'data_25_50_2')
data_folder = os.path.join(monitoring_folder, 'data')
data_8h = os.path.join(data_folder, '8h')

import_path = data_8h

import_dict = {'temp_data': ['arrTemps.txt', ['bt' + str(i) for i in range(1, 20)]],  # change to arrBT
               'mid_data': ['arrMID.txt', ['mid' + str(i) for i in range(1, 7)]],
               'time_data': ['arrTime.txt', ['unix']],
               'q_data': ['arrQ.txt', ['q1', 'q2', 'q3', 'q4', 'q5', 'q6']],  # add 'cop'
               'sw_data': ['arrSW.txt', ['sw4', 'sw14', 'sw17']],
               'sp_data': ['arrSPs.txt', ['mm1', 'mm2', 'mm3', 'mm4', 'mm5', 'gp10', 'gp11', 'ea21']],
               'p_data': ['arrPower.txt', ['p10', 'p11', 'p21', 'p1']]
               }
print(import_dict.keys())


listwithdfs = [pd.read_csv(os.path.join(import_path, value[0]), sep='\t', header=0, names=value[1]) for value in import_dict.values()]

#print(listwithdfs)

temp_data, mid_data, time_data, q_data, sw_data, sp_data, p_data = listwithdfs


'''
#Old import

temp_path = os.path.join(data_8h, 'arrBT.txt')
mid_path = os.path.join(data_8h, 'arrMID.txt')
time_path = os.path.join(data_8h, 'arrTime.txt')
q_path = os.path.join(data_8h, 'arrQ.txt')
sw_path = os.path.join(data_8h, 'arrSW.txt')
sp_path = os.path.join(data_8h, 'arrPID.txt')
p_path = os.path.join(data_8h, 'arrPower.txt')

header_temp = ['bt'+str(i) for i in range(1, 20)]
header_mid = ['mid'+str(i) for i in range(1, 7)]
header_q = ['q'+str(i) for i in range(1, 7)]
header_sw = ['sw4', 'sw14', 'sw17']
header_sp = ['mm1', 'mm2', 'mm3', 'mm4', 'mm5', 'gp10', 'gp11', 'ea21']
header_p = ['p10', 'p11', 'p21', 'p1']

temp_data = pd.read_csv(temp_path, sep='\t', header=0, names=header_temp)
mid_data = pd.read_csv(mid_path, sep='\t', header=0, names=header_mid)
time_data = pd.read_csv(time_path)
time_data['Europe/Berlin'] = pd.to_datetime(time_data['# Unix Time'], unit='s').astype('datetime64[ns, Europe/Berlin]')
q_data = pd.read_csv(q_path, sep='\t', header=0, names=header_q)*1000
sw_data = pd.read_csv(sw_path, sep='\t', header=0, names=header_sw)
sp_data = pd.read_csv(sp_path, sep='\t', header=0, names=header_sp)
p_data = pd.read_csv(p_path, sep='\t', header=0, names=header_p)
'''

"""
Test values
"""
ua = 456
sr = 3000  # Sample row number
obs_w = 3600  # Observation window
offset = obs_w / 2

pipe_AK = Pipe(c2k(temp_data.bt5[sr]), 25, 3, 'copper', 'tubolit')
print(pipe_AK.t_in, pipe_AK.t_out, pipe_AK.Q_loss)


# print(get_air_property_at1bar('pr', 20))
# print(get_k_tubolit(60))


class TableResults(EpsNtu):
    def __init__(self, ua, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in):
        """
        This class creates a dictionary and subsequently a pandas data frame with results.

        :param ua:  value of the heat exchanger [W/K]
        :param tc_in_meas: Measured temperature at the colder circuit inlet [k]
        :param tc_out_meas: Measured temperature at the colder circuit outlet [k]
        :param dvc_in: Measured volumetric flow rate of the colder fluid at inlet [m3/s]
        :param th_in_meas: Measured temperature at the hotter circuit inlet [K]
        :param th_out_meas: Measured temperature at the hotter circuit outlet [K]
        :param dvh_in: Measured volumetric flow rate of the hotter fluid at inlet [m3/s]
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
                             'Measured heat transfer': [self.Qc_calc, self.Qh_calc],
                             'Imported heat transfer': [q_data.q1[sr] * 1000, None],
                             'Simulated heat transfer': [self.Q_sim, self.Q_sim]}
        self.columns = list(self.results_dict.keys())
        self.index = ['Cold side', 'Hot side']
        self.results_df = pd.DataFrame(self.results_dict, columns=self.columns, index=self.index)


"""
Sandbox
"""

# emul_df = TableResults(ua, temp_data.bt9[sr], temp_data.bt17[sr], mid_data.mid1[sr],
# temp_data.bt12[sr], temp_data.bt13[sr], mid_data.mid2[sr])
# print(emul_df.results_df)
# serv_df = TableResults(3250, t3, t4, mid6, t1, t2, mid5)
# print(serv_df.results_df)
# recooler_df = TableResults(3500, bt15, bt16, V_BF4, bt13, bt14, mid2)
# print(recooler_df.results_df)


'''
Plotting the temperatures
'''
# plt.style.use('seaborn-whitegrid')
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)

temp_data['timestamp'] = np.linspace(1, len(temp_data.bt1), len(temp_data.bt1))
fig = plt.figure(1)
ax = fig.gca()
x = temp_data['timestamp']
y1 = temp_data['bt10']
y2 = temp_data['bt11']
y3 = temp_data['bt12']
y4 = temp_data['bt13']

plt.plot(x, y1, label='bt10')
plt.plot(x, y2, label='bt11')
plt.plot(x, y3, label='bt12')
plt.plot(x, y4, label='bt13')

plt.xlim(0, len(q_data.q1))
plt.xlabel('Time [s]')
plt.ylabel('Temp [°C]')
# plt.title('KK temperature fluctuation comparison')
plt.grid('on')
plt.legend(loc='best')
# plt.xlim(sr-600, sr+600)

# plt.savefig(simulation_folder+'/Fluctuation_bt9-17vsbt10-11.png', dpi=600, bbox_inches='tight',
# facecolor=fig.get_facecolor(), edgecolor='none')
# plt.show()

'''
Plotting the control parameters
'''
fig = plt.figure(2)
ax = fig.gca()
x = temp_data['timestamp']
y1 = temp_data['bt4']
y2 = sw_data['sw4']
y3 = sp_data['mm1']

plt.plot(x, y1, label='bt4')
plt.plot(x, y2, label='sw4')
plt.plot(x, y3, label='mm1')

plt.xlim(0, len(temp_data.bt14))
plt.xlabel('Time [s]')
plt.ylabel('Temp [°C]/ Signal [%]')
plt.title('AK regulation')
plt.grid('on')
plt.legend(loc='best')
plt.xlim(sr - offset, sr + offset)

# plt.savefig(simulation_folder+'/AK_regulation_old', dpi=600, bbox_inches='tight',
#            facecolor=fig.get_facecolor(), edgecolor='none')
plt.show()

'''
Plotting the heat transfer
'''
qsim_list = [round(EpsNtu(ua, temp_data.bt10[i], temp_data.bt11[i], mid_data.mid1[i], temp_data.bt12[i],
                          temp_data.bt13[i], mid_data.mid2[i]).Q_sim) for i in range(0, len(temp_data.bt1))]
qcalc_list = [round(EpsNtu(ua, temp_data.bt10[i], temp_data.bt11[i], mid_data.mid1[i], temp_data.bt12[i],
                           temp_data.bt13[i], mid_data.mid2[i]).Qc_calc) for i in range(0, len(temp_data.bt1))]

q_data['timestamp'] = np.linspace(1, len(q_data.q1), len(q_data.q1))
q_data['kk_sim'] = qsim_list
q_data['kk_calc'] = qcalc_list
rmse = get_rmse(q_data.kk_sim, q_data.q1)

fig = plt.figure(3)
# plt.plot(q_data.timestamp, q_data.q1, label='imported', linewidth=0.5)
plt.plot(q_data.timestamp, q_data.kk_sim, label='eps-ntu', linewidth=0.5)
plt.plot(q_data.timestamp, q_data.kk_calc, label='calculated', linewidth=0.5)

# plt.xlim(sr-600, sr+600)
plt.xlim(0, len(q_data.q1))
plt.ylim(1000, 15000)
plt.xlabel('Time [s]')
plt.ylabel('Qdot [W]')
plt.title(f'Qdot over bt9-17 (UA = {ua}, RMSE = {rmse})')
plt.grid('on')
plt.legend(loc='best')
# fig.savefig(simulation_folder+f'/Qdot_bt9-17_ua{ua}_rmse{rmse}'+'.png', dpi=600)
# fig.savefig(simulation_folder+f'/Qdot_bt9-17_ua{ua}_rmse{rmse}_zoomed'+'.png', dpi=600)
# plt.show()
