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

        self.bg_used = False
        if bg_dict:
            self.bg_used = True

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
            self.draw_data_expectation(data_hists[index], total_expectation_hists[index], x_axis_label)

    def draw_data_expectation(self, data, expectation, x_axis_label=''):
        self.plotter.add_hist(data, color='black', xerr=True, histtype='errorbar')
        self.plotter.add_hist(expectation, is_denominator=True)

        self.plotter.create_subplots(2, 1, figsize=(10, 10), height_ratios=[1, 0.3])
        self.plotter.set_current_axis(0, 0)
        self.plotter.draw()
        self.plotter.set_yaxis_config(log_scale=True, set_min=1e-1)
        self.plotter.set_xaxis_config(remove_tick_labels=True)
        self.plotter.set_axis_config()
        self.plotter.set_experiment_label()

        self.plotter.set_current_axis(1, 0)
        self.plotter.draw_ratio()
        self.plotter.set_yaxis_config(set_min=0.5, set_max=1.5)
        self.plotter.get_axis(1, 0).set_xlabel(x_axis_label, fontsize=30)
        self.plotter.save_fig()

    def draw(self):
        self.plotter.create_subplots(1, 1, figsize=(10, 10))
        self.plotter.set_current_axis(0, 0)
        self.plotter.draw()
        self.plotter.set_yaxis_config(label='Acc. x Eff.', set_min=0.0, set_max=1.05, set_grid=True)
        self.plotter.set_xaxis_config(set_grid=True)
        self.plotter.set_axis_config()
        self.plotter.set_experiment_label('Simulation')
        self.plotter.save_fig()

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

        self.unfolders[channel + period] = ISRUnFolders(channel, period,
                                                        self.mass_edges,
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

    def draw_acceptance_x_efficiency_plots(self, channel, period, draw_mass=True):
        # acceptance x efficiency
        # from matrix, get unfolded hist and divide by full_phase hist
        channel, period, hist_path_postfix = self.unfolders[channel + period].get_unfolders_config()
        self.plotter.set_period_channel(period, channel, hist_path_postfix)

        unfolded_isr_hist = self.unfolders[channel + period].get_mc_isr_hists('unfolded', self.hist_producer)
        full_phase_isr_hist = self.unfolders[channel + period].get_mc_isr_hists('full_phase', self.hist_producer)

        if draw_mass:
            unfolded = unfolded_isr_hist.get_mass_hist()
            full_phase = full_phase_isr_hist.get_mass_hist()

            hist = unfolded.divide(full_phase)
            x_axis_label = '$\mathit{m ^{' + channel + '}}$ [GeV]'

            self.plotter.add_hist(hist, color='black', xerr=True, histtype='errorbar')
            self.draw()
        else:
            unfolded = unfolded_isr_hist.get_pt_hist(-1)
            full_phase = full_phase_isr_hist.get_pt_hist(-1)

            x_axis_label = '$\mathit{p_T ^{' + channel + '}}$ [GeV]'
            for index in range(len(unfolded)):
                hist = unfolded[index].divide(full_phase[index])

                self.plotter.add_hist(hist, color='black', xerr=True, histtype='errorbar')
                self.draw()

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

    def draw_unfolded_plots(self, channel, period, show_mc=True, bin_width_norm=True, draw_mass=True):
        period = self.unfolders[channel + period].period
        channel = self.unfolders[channel + period].channel
        hist_path_postfix = self.unfolders[channel + period].hist_path_postfix
        self.plotter.set_period_channel(period, channel, hist_path_postfix)

        mc_isr_hist = self.unfolders[channel + period].get_mc_isr_hists('unfolded', self.hist_producer)
        data_isr_hist = self.unfolders[channel + period].get_isr_hists('unfolded')

        if draw_mass:
            unfolded_data = data_isr_hist.get_mass_hist(bin_width_norm)
            expectation = mc_isr_hist.get_mass_hist(bin_width_norm)

            self.draw_data_expectation(unfolded_data, expectation)
        else:

            unfolded_data = data_isr_hist.get_pt_hist(-1, bin_width_norm)
            expectation = mc_isr_hist.get_pt_hist(-1, bin_width_norm)

            x_axis_label = '$\mathit{p_T ^{' + channel + '}}$ [GeV]'

            for index in range(len(unfolded_data)):
                self.draw_data_expectation(unfolded_data[index], expectation[index])

    def draw_uncertainty_plot(self, channel, period):
        hist_path_postfix = self.unfolders[channel + period].hist_path_postfix
        isr_hists = self.unfolders[channel + period].get_isr_hists("full_phase")
        mass_df, pt_df = isr_hists.get_isr_dataframe()

        self.plotter.set_period_channel(period, channel, hist_path_postfix)
        self.plotter.draw_isr_error_df(pt_df)
