from thermodyn.functions import c2k, lpm2cmps, get_cp, get_rho


class Node:
    def __init__(self, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in):
        """
        This class calculates the heat transfer via measured values.

        :param tc_in_meas: Measured temperature at the colder circuit inlet [k]
        :param tc_out_meas: Measured temperature at the colder circuit outlet [k]
        :param dvc_in: Measured volumetric flow rate of the colder fluid at inlet [m3/s]
        :param th_in_meas: Measured temperature at the hotter circuit inlet [K]
        :param th_out_meas: Measured temperature at the hotter circuit outlet [K]
        :param dvh_in: Measured volumetric flow rate of the hotter fluid at inlet [m3/s]
        """
        self.tc_in_meas = c2k(tc_in_meas)
        self.tc_out_meas = c2k(tc_out_meas)
        self.dtc_meas = self.tc_out_meas - self.tc_in_meas
        self.dVc_in = lpm2cmps(dvc_in)
        self.cpc_in = get_cp(self.tc_in_meas)
        self.rhoc_in = get_rho(self.tc_in_meas)
        self.dmc_in = self.rhoc_in * self.dVc_in
        self.Qc_calc = round(self.dmc_in * self.cpc_in * self.dtc_meas)

        self.th_in_meas = c2k(th_in_meas)
        self.th_out_meas = c2k(th_out_meas)
        self.dth_meas = self.th_in_meas - self.th_out_meas
        self.dVh_in = lpm2cmps(dvh_in)
        self.cph_in = get_cp(self.th_in_meas)
        self.rhoh_in = get_rho(self.th_in_meas)
        self.dmh_in = self.rhoh_in * self.dVh_in
        self.Qh_calc = round(self.dmh_in * self.cph_in * self.dth_meas)