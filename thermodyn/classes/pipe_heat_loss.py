from thermodyn.functions import get_cp, get_rho, get_k, c2k, get_t_avg, get_surf_area_cyl, get_cond_res

"""
Duct sizes dictionary for heat loss calculation
"""
dict_pipe_sizes = {25: {'ri': 0.0136, 'ro': 0.01685},  # DIN EN 10255
                   32: {'ri': 0.01795, 'ro': 0.0212},
                   40: {'ri': 0.0209, 'ro': 0.02415}}


class Pipe:
    def __init__(self, t_in, dn, length_pipe, mat_pipe, mat_insul, ho=7, t_amb=22):
        """
        This class calculates the output temperature of the pipe.

        :param t_in: inlet temperature [K]
        :param t_amb: ambient temperature [K]
        :param r_pi: inner radius of the pipe [m]
        :param r_po: outer radius of the pipe [m]
        :param r_io: outer radius of the insulated pipe [m]
        :param h_o: heat transfer coefficient between the outer surface of the insulated pipe and ambient [W/m2K]
        #https://local.armacell.com/fileadmin/cms/uk/service/en/ArmaflexUKSpecGuideDigital.pdf
        :param k_p: thermal conductivity of the pipe material [W/mK]
        :param k_i: thermal conductivity of the insulation material [W/mK]
        """
        self.t_amb = c2k(t_amb)
        self.L = length_pipe
        self.ho = ho
        self.t_in = t_in
        self.t_air_film = (self.t_in + self.t_amb)/2
        self.dn = dn
        self.r_pi = dict_pipe_sizes[dn]['ri']
        self.r_po = dict_pipe_sizes[dn]['ro']
        self.r_io = self.r_po + 0.009
        self.k_p = get_k(mat_pipe, self.t_in)
        self.k_i = get_k(mat_insul, self.t_air_film)

        self.A1 = get_surf_area_cyl(self.r_pi, self.L)
        self.A3 = get_surf_area_cyl(self.r_io, self.L)

        # self.Ri = 1 / (self.h1 * self.A1) negligibly small?
        self.R1 = get_cond_res(self.r_po, self.r_pi, self.k_p, self.L)
        self.R2 = get_cond_res(self.r_io, self.r_po, self.k_i, self.L)
        self.Ro = 1 / (self.ho * self.A3)
        self.R_tot = self.R1 + self.R2 + self.Ro

        #self.t_out = self.t_in
        #for dL in np.linspace(0, self.L, 10):
        #    self.Q_loss = (self.t_out - self.t_amb) / self.R_tot
        #    self.t_out = self.t_in - self.Q_loss * dL / (get_cp(self.t_out) * get_rho(self.t_out) * self.A1 * dL)

        self.Q_loss = (self.t_in - self.t_amb) / self.R_tot
        self.t_out = self.t_in - self.Q_loss * self.L / (get_cp(self.t_in) * get_rho(self.t_in) * self.A1 * self.L)

