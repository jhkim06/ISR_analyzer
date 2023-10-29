import ctypes

import ROOT
from collections import namedtuple
from types import MappingProxyType
from Hist import Hist
import re
import numpy as np


State = namedtuple('State', ['RegMode', 'DensityMode', 'EConstraint'])


class UnFolder:
    def __init__(self, response_matrix, input_hist, bg_hists=None,
                 reg_mode="None", ex_constraint="Area", density_mode="None"):
        self.reg_mode = MappingProxyType(
            {"None": ROOT.TUnfold.kRegModeNone,
             "Size": ROOT.TUnfold.kRegModeSize,
             "Derivative": ROOT.TUnfold.kRegModeDerivative,
             "Curvature": ROOT.TUnfold.kRegModeCurvature,
             "Mixed": ROOT.TUnfold.kRegModeMixed
             }
        )
        # TODO study custom regularisation
        self.density_mode = MappingProxyType(
            {"None": ROOT.TUnfoldDensity.kDensityModeNone,
             "BinWidth": ROOT.TUnfoldDensity.kDensityModeBinWidth,
             "User": ROOT.TUnfoldDensity.kDensityModeUser,
             "BinWidthAndUser": ROOT.TUnfoldDensity.kDensityModeBinWidthAndUser
             }
        )
        self.ex_constraint = MappingProxyType(
            {"None": ROOT.TUnfold.kEConstraintNone,
             "Area": ROOT.TUnfold.kEConstraintArea
             }
        )

        self.response_matrix = response_matrix
        # 2D unfolding use TUnfoldBinning
        self.folded_bin = None
        self.unfolded_bin = None
        self.tunfold_bin_used = False
        self.use_axis_binning = False
        self.set_bins()

        # input hist
        self.input_hist = input_hist
        # background hists
        self.bg_hists = bg_hists

        self.default_tunfold_density = self.create_tunfold_density(reg_mode, ex_constraint, density_mode)
        self.set_input(self.default_tunfold_density)  # set input histogram, and subtract backgrounds
        # for systematics
        self.sys_tunfold_density = dict()
        self.has_systematic = False
        self.set_systematics(reg_mode, ex_constraint, density_mode)

        self.reg_reg_strength = 0

    def get_projected_hist(self, axis='x', first_bin=0, last_bin=-1, option='e'):
        return self.response_matrix.get_projected_hist(axis, first_bin, last_bin, option)

    def set_bins(self):
        # get folded bin from the input_hist
        if "[tunfold-matrix]" in self.response_matrix.hist_name:
            folded_bin, unfolded_bin = self.response_matrix.get_bin_name()
            self.folded_bin = self.response_matrix.get_unfold_bin(folded_bin)
            self.unfolded_bin = self.response_matrix.get_unfold_bin(unfolded_bin)
            self.tunfold_bin_used = True
        else:
            self.tunfold_bin_used = False
            self.use_axis_binning = True

    def set_systematic(self, reg_mode, ex_constraint, density_mode, root_sys_hist_dict):
        self.has_systematic = True
        for sys_name, variations in root_sys_hist_dict.items():
            if sys_name in self.sys_tunfold_density.keys():
                continue
            self.sys_tunfold_density[sys_name] = dict()
            for variation in variations:
                self.sys_tunfold_density[sys_name][variation] = \
                    self.create_tunfold_density(reg_mode, ex_constraint, density_mode, sys_name, variation)
                self.set_input(self.sys_tunfold_density[sys_name][variation], sys_name, variation)

    def set_systematics(self, reg_mode, ex_constraint, density_mode):

        if self.input_hist.root_sys_hist_dict is not None:
            self.set_systematic(reg_mode, ex_constraint, density_mode, self.input_hist.root_sys_hist_dict)
        if self.response_matrix.root_sys_hist_dict is not None:
            self.set_systematic(reg_mode, ex_constraint, density_mode, self.response_matrix.root_sys_hist_dict)
        if self.bg_hists is not None and self.bg_hists.root_sys_hist_dict is not None:
            self.set_systematic(reg_mode, ex_constraint, density_mode, self.bg_hists.root_sys_hist_dict)

    def create_tunfold_density(self, reg_mode, ex_constraint, density_mode,
                               sys_name="", sys_variation=""):

        if self.tunfold_bin_used:
            tunfold_density = ROOT.TUnfoldDensity(self.response_matrix.get_root_hist(sys_name, sys_variation),
                                                  ROOT.TUnfold.kHistMapOutputHoriz,
                                                  self.reg_mode[reg_mode],
                                                  self.ex_constraint[ex_constraint],
                                                  self.density_mode[density_mode],
                                                  self.unfolded_bin, self.folded_bin)
        else:
            # print("use 1D hists...")
            tunfold_density = ROOT.TUnfoldDensity(self.response_matrix.get_root_hist(sys_name, sys_variation),
                                                  ROOT.TUnfold.kHistMapOutputHoriz,
                                                  self.reg_mode[reg_mode],
                                                  self.ex_constraint[ex_constraint],
                                                  self.density_mode[density_mode])

        return tunfold_density

    def set_input(self, tunfold_density, sys_name="", sys_variation=""):
        # set input histogram and subtract background histograms
        tunfold_density.SetInput(self.input_hist.get_root_hist(sys_name, sys_variation))
        fake_hist_name = re.sub(r'reco', 'reco_fake', self.input_hist.hist_name)
        fake_hist = self.response_matrix.file_group.get_hist(fake_hist_name)
        tunfold_density.SubtractBackground(fake_hist.get_root_hist(sys_name, sys_variation), "DY fake")

        if self.bg_hists is not None:
            for label, hist in self.bg_hists.hist_dict.items():
                tunfold_density.SubtractBackground(hist.get_root_hist(sys_name, sys_variation), label)

    def do_unfold(self, unfold_method=""):
        # TODO record unfold state and use the same state for systematic unfolding
        if unfold_method == "":  # chi-square minimization
            self.reg_reg_strength = 0
            self.default_tunfold_density.DoUnfold(self.reg_reg_strength)
        elif unfold_method == "scan_tau":
            self.reg_reg_strength = 0  # TODO update strength and save it
        elif unfold_method == "scan_lcurve":
            self.reg_reg_strength = 0
        elif unfold_method == "partial_reg":
            self.reg_reg_strength = 0
        else:
            return
        self.do_sys_unfolds()
        # self.bottom_line_test()

    def do_sys_unfolds(self):
        for sys_name, variations in self.sys_tunfold_density.items():
            for variation in variations:
                self.sys_tunfold_density[sys_name][variation].DoUnfold(self.reg_reg_strength)

    def get_unfolded_hist_name(self):
        res_name = self.response_matrix.hist_name
        matches = re.findall(r'(\[[^]]*])', self.response_matrix.hist_name)
        if self.tunfold_bin_used:
            hist_name = "[tunfold-hist]_" + matches[1] + "_" + matches[3]
        else:
            hist_name = res_name.split(matches[0])[0] + matches[1] + res_name.split(matches[1])[-1]
        return hist_name

    def get_folded_hist_name(self):
        res_name = self.response_matrix.hist_name
        matches = re.findall(r'(\[[^]]*])', self.response_matrix.hist_name)
        if self.tunfold_bin_used:
            hist_name = "[tunfold-hist]_" + matches[1] + "_" + matches[2]
        else:
            hist_name = res_name.split(matches[0])[0] + matches[0] + res_name.split(matches[1])[-1]
        return hist_name

    def get_unfolded_hist(self):
        # change histogram name to contain the unfolded bin name
        hist_name = self.get_unfolded_hist_name()
        unfolded_hist = self.get_tunfold_output()

        sys_hist_dict = None
        if self.has_systematic:
            sys_hist_dict = dict()
            for sys_name, variations in self.sys_tunfold_density.items():
                sys_hist_dict[sys_name] = dict()
                for variation in variations:
                    sys_hist_dict[sys_name][variation] = self.get_tunfold_output(sys_name, variation)

        hist = Hist(hist_name, unfolded_hist, sys_hist_dict, self.response_matrix.file_group)
        hist.is_mc = self.input_hist.is_mc
        hist.update_group_label("Unfolded Data")
        return hist

    def get_tunfold_output(self, sys_name="", variation=""):
        if sys_name == "":
            return self.default_tunfold_density.GetOutput("default",
                                                          ctypes.c_char_p(0),
                                                          ctypes.c_char_p(0), "*[*]",
                                                          self.use_axis_binning)
        else:
            return self.sys_tunfold_density[sys_name][variation].GetOutput(sys_name + "_" + variation,
                                                                           ctypes.c_char_p(0),
                                                                           ctypes.c_char_p(0), "*[*]",
                                                                           self.use_axis_binning)

    def get_folded_hist(self):
        hist_name = self.get_folded_hist_name()
        folded_hist = self.get_tunfold_input()

        sys_hist_dict = None
        if self.has_systematic:
            sys_hist_dict = dict()
            for sys_name, variations in self.sys_tunfold_density.items():
                sys_hist_dict[sys_name] = dict()
                for variation in variations:
                    sys_hist_dict[sys_name][variation] = self.get_tunfold_input(sys_name, variation)

        hist = Hist(hist_name, folded_hist, sys_hist_dict, self.input_hist.file_group)
        hist.is_mc = self.input_hist.is_mc
        hist.update_group_label("Unfolded Data")
        return hist

    def get_tunfold_input(self, sys_name="", variation=""):  # return
        if sys_name == "":
            return self.default_tunfold_density.GetInput("unfold_input",
                                                         ctypes.c_char_p(0),
                                                         ctypes.c_char_p(0), "*[*]",
                                                         self.use_axis_binning)
        else:
            return self.sys_tunfold_density[sys_name][variation].GetInput(sys_name + "_" + variation,
                                                                          ctypes.c_char_p(0),
                                                                          ctypes.c_char_p(0), "*[*]",
                                                                          self.use_axis_binning)

    def get_chi_square(self, is_unfolded=True):
        if is_unfolded:
            data_hist = self.get_unfolded_hist()
            mc_hist = self.get_expectation_hist(unfolded=True, use_matrix=True)
            text = 'unfolded'
        else:
            data_hist = self.get_folded_hist()  # Todo check if this is bg subtracted data hist
            mc_hist = self.get_expectation_hist(unfolded=False, use_matrix=True)
            text = 'folded'

        data_hist.set_hist_config(True, '')  # note that bin_width_norm forced to be True
        mc_hist.set_hist_config(True, '')
        chi2 = np.sum(np.square(data_hist.to_numpy().values - mc_hist.to_numpy().values) /
                      mc_hist.to_numpy().values)
        print("---bottom_line_test---")
        print(f'{text} {chi2}')

        if self.tunfold_bin_used:
            data_hists = data_hist.get_1d_hists()
            mc_hists = mc_hist.get_1d_hists()
            # for each window
            for index in range(len(data_hists)):
                data_hist_ = data_hists[index]  # Todo check if this is bg subtracted data hist
                mc_hist_ = mc_hists[index]
                chi2 = np.sum(np.square(data_hist_.to_numpy().values - mc_hist_.to_numpy().values) /
                              mc_hist_.to_numpy().values)
                print(f'{text} {index} {chi2}', end=' ')
        print()

    def bottom_line_test(self):

        self.get_chi_square(is_unfolded=False)
        self.get_chi_square()

    # get expectation histogram
    def get_expectation_hist(self, unfolded=True, use_matrix=False):
        if use_matrix:  # projection of response matrix
            if unfolded:
                unfolded_hist = self.get_projected_hist(axis='x')
                return unfolded_hist
            else:
                folded_hist = self.get_projected_hist(axis='y')
                return folded_hist
        else:  # get histogram from the root file
            if unfolded:
                hist_name = self.get_unfolded_hist_name()
                unfolded_hist = self.response_matrix.file_group.get_hist(hist_name)
                return unfolded_hist
            else:
                # folded mc
                hist_name = self.get_folded_hist_name()
                folded_hist = self.response_matrix.file_group.get_hist(hist_name)
                return folded_hist

    def scan_tau(self):
        pass

    def scan_lcurve(self):
        pass
