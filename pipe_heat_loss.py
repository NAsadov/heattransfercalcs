from math import pi


"""
Duct sizes dictionary for heat loss calculation
"""
dict_pipe_sizes = {25: {'ri': 0.0136, 'ro': 0.01685},  # DIN EN 10255
                   32: {'ri': 0.01795, 'ro': 0.0212},
                   40: {'ri': 0.0209, 'ro': 0.02415}}


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
