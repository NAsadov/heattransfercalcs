from CoolProp.CoolProp import PropsSI
from functions import c2k, lpm2cmps

'''
class ThermodynamicProperties:
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
        self.tc_in_meas = c2k(tc_in_meas)
        self.tc_out_meas = c2k(tc_out_meas)
        self.dtc_meas = self.tc_out_meas - self.tc_in_meas
        self.dVc_in = lpm2cmps(dvc_in)
        self.cpc_in = PropsSI('C', 'T', self.tc_in_meas, 'P', p, fluid)
        self.rhoc_in = PropsSI('D', 'T', self.tc_in_meas, 'P', p, fluid)
        self.dmc_in = self.rhoc_in * self.dVc_in
        self.Qc_calc = round(self.dmc_in * self.cpc_in * self.dtc_meas)

        self.th_in_meas = c2k(th_in_meas)
        self.th_out_meas = c2k(th_out_meas)
        self.dth_meas = self.th_in_meas - self.th_out_meas
        self.dVh_in = lpm2cmps(dvh_in)
        self.cph_in = PropsSI('C', 'T', self.th_in_meas, 'P', p, fluid)
        self.rhoh_in = PropsSI('D', 'T', self.th_in_meas, 'P', p, fluid)
        self.dmh_in = self.rhoh_in * self.dVh_in
        self.Qh_calc = round(self.dmh_in * self.cph_in * self.dth_meas)