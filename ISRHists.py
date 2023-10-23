import pandas as pd
import numpy as np


def calculate_squared_root_sum(raw_df, reg_expression, new_col_name="total error"):

    selected_cols = raw_df.filter(regex=reg_expression)
    squared_sum = (selected_cols ** 2).sum(axis=1)
    root_of_squared_sum = np.sqrt(squared_sum)

    raw_df[new_col_name] = root_of_squared_sum


class ISRHists:
    def __init__(self, mass_hist, *pt_hist, mass_windows=None):  # mass windows

        self.tunfold_bins = dict()
        self.mass_windows = mass_windows

        # TODO use ISRHists as input to unfold
        self.mass_2d_hist = None
        self.pt_2d_hist = None

        # always one mass histogram
        if "tunfold" in mass_hist.hist_name:
            self.mass_2d_hist = mass_hist
            self.mass_isr_hist = mass_hist.get_1d_hists()[0]
        else:
            self.mass_isr_hist = mass_hist

        # pt
        if len(pt_hist) == 1:  # FIXME use explicit flag for 2D hist
            self.pt_2d_hist = pt_hist[0]
            self.pt_isr_hists = self.pt_2d_hist.get_1d_hists()

            bin_name, _ = self.pt_2d_hist.get_bin_name()
            pt_unfold_bin = self.pt_2d_hist.get_unfold_bin(bin_name)
            edges = pt_unfold_bin.GetDistributionBinning(1)
            self.mass_windows = [(edges[i], edges[i+1]) for i in range(len(edges)-1)]
        else:
            self.pt_isr_hists = pt_hist

    def get_pt_hist(self, mass_index, bin_width_norm=False):
        # TODO error handle for invalid mass_index
        if mass_index == -1:
            for index, hist in enumerate(self.pt_isr_hists):
                hist.set_hist_config(bin_width_norm, '')
                text = ("{:.0f}".format(float(self.mass_windows[index][0])) + r"$\mathit{ < m < }$" +
                        "{:.0f}".format(float(self.mass_windows[index][1])) + " GeV")
                hist.set_text_for_plot(text)
            return self.pt_isr_hists
        else:
            self.pt_isr_hists[mass_index].set_hist_config(bin_width_norm, '')
            return self.pt_isr_hists[mass_index]

    def get_mass_hist(self, bin_width_norm=False):
        self.mass_isr_hist.set_hist_config(bin_width_norm, '')
        return self.mass_isr_hist

    def get_isr_dataframe(self, binned_mean=False):
        # get pt, mass mean for each mass window
        dipt_dict_list = []
        dimass_dict_list = []

        for ith_mass in range(len(self.mass_windows)):
            # get mean mass and pt
            low_mass_edge = self.mass_windows[ith_mass][0]
            high_mass_edge = self.mass_windows[ith_mass][1]
            mean_mass, mean_mass_error = self.mass_isr_hist.get_mean(low_mass_edge, high_mass_edge)
            mean_pt, mean_pt_error = self.pt_isr_hists[ith_mass].get_mean(binned_mean=binned_mean)

            mass_window_str = str(low_mass_edge) + ":" + str(high_mass_edge)
            dipt_dict = {
                "mass_window": mass_window_str,
                "mean": mean_pt}
            dimass_dict = {
                "mass_window": mass_window_str,
                "mean": mean_mass
            }

            for sys_name in mean_mass_error.keys():
                dimass_dict[sys_name + " error"] = mean_mass_error[sys_name]
            for sys_name in mean_pt_error.keys():
                dipt_dict[sys_name + " error"] = mean_pt_error[sys_name]
            # error
            dipt_dict_list.append(dipt_dict)
            dimass_dict_list.append(dimass_dict)

        dimass_df = pd.DataFrame(dimass_dict_list)
        dipt_df = pd.DataFrame(dipt_dict_list)

        # calculate total error
        reg_expression = r".*error$"
        calculate_squared_root_sum(dimass_df, reg_expression)
        calculate_squared_root_sum(dipt_df, reg_expression)

        return dimass_df, dipt_df
