from math import exp
from basic_calc import Node


class EpsNtu(Node):
    def __init__(self, ua, tc_in_meas, tc_out_meas, dvc_in, th_in_meas, th_out_meas, dvh_in):
        """
        This class calculates the heat transfer and the output temperatures via epsilon-NTU method.

        :param ua:  value of the heat exchanger [W/K]
        :param tc_in_meas: Measured temperature at the colder circuit inlet [k]
        :param tc_out_meas: Measured temperature at the colder circuit outlet [k]
        :param dvc_in: Measured volumetric flow rate of the colder fluid at inlet [m3/s]
        :param th_in_meas: Measured temperature at the hotter circuit inlet [K]
        :param th_out_meas: Measured temperature at the hotter circuit outlet [K]
        :param dvh_in: Measured volumetric flow rate of the hotter fluid at inlet [m3/s]
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
