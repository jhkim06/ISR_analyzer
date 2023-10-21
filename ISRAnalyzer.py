from Analyzer import Analyzer
from ISRHists import ISRHists
import re
from ISRUnFolders import ISRUnFolders


def get_full_phase_hist_name(hist):
    hist_name = re.sub(r'(__)', '_acceptance__', hist.hist_name)
    return hist_name


class ISRAnalyzer(Analyzer):
    def __init__(self, file_pather, signal_dict=None, bg_dict=None):

        super().__init__(file_pather, signal_dict, bg_dict)
        self.analysis_name = "ISR"

        # major histogram name for ISR analysis
        self.pt_hist_name = ''
        self.mass_hist_name = ''
        self.pt_matrix_name = ''
        self.mass_matrix_name = ''

        self.mass_edges = None  # pair of bin edges for 1D setting, used to make mass_postfix_for_pt_hist
        self.pt_edges = None
        self.mass_postfix_for_pt_hist = None
        self.pt_postfix_for_mass_hist = None
        self.is_2d = True  # if 2D, dipt histograms for each mass window in one 1D histogram

        self.unfolders = {}

        self.bg_used = True
        if bg_dict is None:
            self.bg_used = False

        self.input_data_isr_hists = None  # before bkg. subtraction

    def set_mass_pt_hist_names(self, mass_hist_name, pt_hist_name,
                               mass_matrix_name, pt_matrix_name,
                               mass_edges=None, pt_edges=None):

        self.mass_hist_name = mass_hist_name
        self.mass_matrix_name = mass_matrix_name
        self.pt_hist_name = pt_hist_name
        self.pt_matrix_name = pt_matrix_name

        if mass_edges is not None and pt_edges is not None:
            self.is_2d = False
            self.mass_edges = mass_edges
            self.pt_edges = pt_edges

            self.mass_postfix_for_pt_hist = [
                '_dimass_{}to{}'.format(float(mass_edges[index][0]), float(mass_edges[index][1]))
                for index in range(len(mass_edges))
            ]
            self.pt_postfix_for_mass_hist = [
                '_dipt_{}to{}'.format(float(pt_edges[index][0]), float(pt_edges[index][1]))
                for index in range(len(pt_edges))
            ]
        else:
            self.is_2d = True

    def draw_detector_plots(self, channel, period, file_path_postfix="", hist_path_postfix="",
                            draw_mass=True, bin_width_norm=True):
        self.hist_producer.set_base_configs(period, channel, file_path_postfix, hist_path_postfix)
        if draw_mass:
            hist_name = self.mass_hist_name
            second_axis_postfix = self.pt_postfix_for_mass_hist
            x_axis_label = '$\mathit{m ^{' + channel + '}}$ [GeV]'
        else:
            hist_name = self.pt_hist_name
            second_axis_postfix = self.mass_postfix_for_pt_hist
            x_axis_label = '$\mathit{p_T ^{' + channel + '}}$ [GeV]'
        self.hist_producer.set_configs(hist_name,
                                       bin_width_norm=bin_width_norm,
                                       second_axis_postfix=second_axis_postfix,
                                       set_text=True)
        if self.is_2d:
            force_1d_output = True
        else:
            force_1d_output = False

        data_hists = self.hist_producer.get_data_hist(force_1d_output=force_1d_output)
        total_expectation_hists = self.hist_producer.get_total_expectation_hist(force_1d_output=force_1d_output)

        self.plotter.set_period_channel(period, channel, hist_path_postfix)
        for index in range(len(data_hists)):
            self.plotter.draw_comparison(total_expectation_hists[index], data_hists[index],
                                         x_axis_label=x_axis_label,
                                         logy=True, ymin=0.5, ymax=1.5)

    def get_unfold_input_hists_and_matrix(self, is_mass=True):
        if is_mass:
            hist_name = self.mass_hist_name
            matrix_name = self.mass_matrix_name
            second_axis_edges = self.pt_postfix_for_mass_hist
        else:
            hist_name = self.pt_hist_name
            matrix_name = self.pt_matrix_name
            second_axis_edges = self.mass_postfix_for_pt_hist

        # set histogram config
        self.hist_producer.set_configs(hist_name,
                                       second_axis_postfix=second_axis_edges)
        # set data histogram
        input_hist = self.hist_producer.get_data_hist()  # list of 1D or one 2D
        # set background histograms
        bg_hist = None
        if self.bg_used:
            bg_hist = self.hist_producer.get_total_expectation_hist(exp_type='bg')
        # set matrix
        self.hist_producer.set_configs(matrix_name,
                                       second_axis_postfix=second_axis_edges)
        matrix = self.hist_producer.get_total_expectation_hist(exp_type='signal')
        return input_hist, bg_hist, matrix

    def set_unfold(self, channel, period, file_path_postfix="", hist_path_postfix=""):
        # TODO allow mass bin option
        # initialize the list of mass and pt unfolder before append
        self.hist_producer.set_base_configs(period, channel, file_path_postfix, hist_path_postfix)
        mass_data_hist, mass_bg_hist, mass_matrix = self.get_unfold_input_hists_and_matrix()
        pt_data_hist, pt_bg_hist, pt_matrix = self.get_unfold_input_hists_and_matrix(is_mass=False)

        self.unfolders[channel + period] = ISRUnFolders(channel, period, self.mass_edges,
                                                        file_path_postfix, hist_path_postfix,
                                                        self.is_2d)
        self.unfolders[channel + period].set_unfolders(mass_data_hist, mass_matrix, mass_bg_hist)
        self.unfolders[channel + period].set_unfolders(pt_data_hist, pt_matrix, pt_bg_hist, unfolder_name='pt')

        if self.is_2d:
            self.input_data_isr_hists = ISRHists(mass_data_hist, pt_data_hist,
                                                 mass_windows=self.mass_edges)
        else:
            self.input_data_isr_hists = ISRHists(mass_data_hist[0], *pt_data_hist,
                                                 mass_windows=self.mass_edges)

    def unfold(self, channel, period):
        self.unfolders[channel + period].unfold()

    def acceptance_corrections(self, channel, period):
        self.unfolders[channel + period].acceptance_corrections(self.hist_producer)

    def draw_isr_plots(self, channel, period, level_name, binned_mean=True, show_mc=False):
        period = self.unfolders[channel + period].period
        channel = self.unfolders[channel + period].channel
        hist_path_postfix = self.unfolders[channel + period].hist_path_postfix
        self.plotter.set_period_channel(period, channel, hist_path_postfix)

        # get isr_hists
        isr_hists = self.unfolders[channel + period].get_isr_hists(level_name)
        data_mass_df, data_pt_df = isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        label = level_name + ' data (' + self.unfolders[channel + period].channel + ')'
        self.plotter.draw_isr_data_frame(data_mass_df, data_pt_df,
                                         ymin=14, ymax=30,
                                         new_fig=True,
                                         color='black', label=label)
        if show_mc:
            isr_hists = self.unfolders[channel + period].get_mc_isr_hists(level_name, self.hist_producer)
            data_mass_df, data_pt_df = isr_hists.get_isr_dataframe(binned_mean=binned_mean)
            label = level_name + ' mc (' + self.unfolders[channel + period].channel + ')'
            self.plotter.draw_isr_data_frame(data_mass_df, data_pt_df,
                                             ymin=14, ymax=30,
                                             new_fig=False,
                                             color='gray', label=label)

    # def draw_unfolded_data(self, channel, period, show_mc=True, bin_width_norm=True, draw_mass=True)
    def draw_unfolded_data(self, channel, period, show_mc=True, bin_width_norm=True, draw_mass=True):
        period = self.unfolders[channel + period].period
        channel = self.unfolders[channel + period].channel
        hist_path_postfix = self.unfolders[channel + period].hist_path_postfix
        self.plotter.set_period_channel(period, channel, hist_path_postfix)

        mc_isr_hist = self.unfolders[channel + period].get_mc_isr_hists('unfolded', self.hist_producer)
        data_isr_hist = self.unfolders[channel + period].get_isr_hists('unfolded')

        if draw_mass:
            unfolded_data = data_isr_hist.get_mass_hist(bin_width_norm)
            expectation = mc_isr_hist.get_mass_hist(bin_width_norm)

            x_axis_label = '$\mathit{m ^{' + channel + '}}$ [GeV]'
            self.plotter.draw_comparison(expectation, unfolded_data,
                                         x_axis_label=x_axis_label)
        else:

            unfolded_data = data_isr_hist.get_pt_hist(-1, bin_width_norm)
            expectation = mc_isr_hist.get_pt_hist(-1, bin_width_norm)

            x_axis_label = '$\mathit{p_T ^{' + channel + '}}$ [GeV]'
            for index in range(len(unfolded_data)):
                self.plotter.draw_comparison(expectation[index], unfolded_data[index],
                                             x_axis_label=x_axis_label)

    def draw_uncertainty_plot(self):
        mass_df, pt_df = self.unfolded_isr_hists.get_isr_dataframe()
        self.plotter.draw_isr_error_df(pt_df)
