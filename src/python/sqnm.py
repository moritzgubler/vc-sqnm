import numpy as np
import historylist

class SQNM:
    def __init__(self, ndim, nhist_max, eps_supsp, alpha_min):
        self.ndim = ndim
        self.nhist_max = nhist_max
        self.eps_subsp = eps_supsp
        self.alpha_min = alpha_min
        self.xlist = historylist.HistoryList(self.ndim, self.nhist_max)
        self.flist = historylist.HistoryList(self.ndim, self.nhist_max)
        self.alpha = 0.0
        self.dir_of_descent = np.zeros(ndim)
        self.prev_f_of_x = 0.0
        self.prev_df_dx = np.zeros(ndim)
        self.s_evec = np.zeros((nhist_max, nhist_max))
        self.s_eval = np.zeros(nhist_max)
        self.dr_subsp = np.zeros((ndim, nhist_max))
        self.df_subsp = np.zeros((ndim, nhist_max))
        self.h_evec_subsp = np.zeros((nhist_max, nhist_max))
        self.h_evec = np.zeros((ndim, nhist_max))
        self.h_eval = np.zeros(nhist_max)
        self.res = np.zeros(nhist_max)
        self.res_temp = np.zeros(ndim)
        self.gainratio = 0.0
        self.nhist = 0

    def sqnm_step(self, x, f_of_x, df_dx):
        pass



