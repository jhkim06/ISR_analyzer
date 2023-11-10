import ROOT
import ctypes
import numpy as np


ROOT.gSystem.Load("/Users/jhkim/cms_snu/2.4.0/libBlue.so")

for_value = "%5.2f"
for_uncertainty = for_value
for_weight = "%6.4f"
for_rho = "%5.2f"
for_pull = "%4.2f"
for_chi = "%5.3f"
for_unit = " GeV"


class Combiner:
    def __init__(self,
                 observable_name,
                 estimation_info,
                 use_rho_same_channel=True,
                 ):
        self.rho_same_channel = {
            "stat": 0.0, "IDsf": 0.0
        }
        self.rho_same_period = {
            "stat": 0.0, "IDsf": 1.0
        }

        self.n_est = len(estimation_info)
        self.n_unc = len(estimation_info[0][1])
        self.n_obs = 1
        self.use_rho_same_channel = use_rho_same_channel

        self.observable_name = observable_name
        self.estimation_info = estimation_info
        self.estimation_values = []
        self.unc_names = []
        self.make_estimation_value_list()
        self.make_unc_name_list()

        self.rho_sou = ROOT.TMatrixD(self.n_est, self.n_est)
        self.blue = ROOT.Blue(self.n_est, self.n_unc)
        self.blue.SetFormat(
            for_value,
            for_uncertainty,
            for_weight,
            for_rho,
            for_pull,
            for_chi,
            for_unit
        )

    def make_estimation_value_list(self):
        for estimation in self.estimation_info:
            self.estimation_values.append(estimation[0][1])
            for uncertainty in estimation[1]:
                self.estimation_values.append(uncertainty[1])

    def make_unc_name_list(self):
        for uncertainty in self.estimation_info[0][1]:
            self.unc_names.append(uncertainty[0])

    def set_name_for_estimation(self):
        name_string = "{"
        for estimation in self.estimation_info:
            name_string = name_string + "\"" + estimation[0][0] + "\","
        name_string += "};"

        cpp_string = f"const TString NamEst[{self.n_est}] = "
        cpp_string += name_string
        cpp_string += "const TString *const p_NamEst = &NamEst[0];"
        print(cpp_string)
        ROOT.gInterpreter.Declare(cpp_string)
        self.blue.FillNamEst(ROOT.p_NamEst)

    def set_name_for_uncertainty(self):
        name_string = "{"
        for uncertainty in self.estimation_info[0][1]:
            name_string = name_string + "\"" + uncertainty[0] + "\","
        name_string += "};"

        cpp_string = f"const TString NamUnc[{self.n_unc}] = "
        cpp_string += name_string
        cpp_string += "const TString *const p_NamUnc = &NamUnc[0];"
        print(cpp_string)
        ROOT.gInterpreter.Declare(cpp_string)
        self.blue.FillNamUnc(ROOT.p_NamUnc)

    def set_name_for_observable(self):
        cpp_string = f"const TString NamObs[1] = "
        cpp_string += "{\"" + self.observable_name + "\"};"
        cpp_string += "const TString *const p_NamObs = &NamObs[0];"
        print(cpp_string)
        ROOT.gInterpreter.Declare(cpp_string)
        self.blue.FillNamObs(ROOT.p_NamObs)

    def set_names(self):
        self.set_name_for_estimation()
        self.set_name_for_uncertainty()
        self.set_name_for_observable()

    def fill_estimation(self):
        estimation_index = 0
        for i in range(self.n_est):
            self.blue.FillEst(i,
                              np.array(self.estimation_values[estimation_index:])
                              )
            print(self.estimation_values[estimation_index:])
            estimation_index += self.n_unc + 1

    def fill_correlation(self):
        for i in range(self.n_unc):
            self.fill_matrix(self.unc_names[i])  # FIXME use name
            self.blue.FillCor(i, self.rho_sou)

    def fill_matrix(self, unc_name):

        if self.use_rho_same_channel:
            rho = self.rho_same_channel
        else:
            rho = self.rho_same_period

        for i in range(self.n_est):
            for j in range(self.n_est):
                if i == j:
                    self.rho_sou[i][j] = 1
                else:
                    self.rho_sou[i][j] = rho[unc_name]

    def set_inputs(self):
        self.fill_estimation()
        self.fill_correlation()

    def solve(self):
        self.blue.FixInp()
        self.blue.Solve()

    def print(self):
        self.blue.PrintEst()
        self.blue.PrintResult()

    def get_result_list(self):
        result_matrix = ROOT.TMatrixD(
            self.blue.GetActObs(),
            self.blue.GetActUnc()+1
        )
        self.blue.GetResult(result_matrix)
        result_list = []
        combined_result = result_matrix[0][0]
        result_list.append(combined_result)

        unc_list = []
        for unc_index in range(1, self.blue.GetActUnc()+1):
            unc_list.append(result_matrix[0][unc_index])

        total_unc = np.sqrt(np.sum(np.square(unc_list)))
        result_list.extend(unc_list)
        result_list.append(total_unc)

        return result_list
