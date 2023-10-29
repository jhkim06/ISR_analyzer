from UnFolder import UnFolder
from Analyzer import Analyzer
from ISRHists import ISRHists
import re


def make_axis_label(channel, is_mass=True):
    if channel == "ee":
        channel_label = "ee"
    elif channel == "mm":
        channel_label = "\mu\mu"
    else:
        channel_label = "ll"

    if is_mass:
        return r'$m^{' + channel_label + '}$ [GeV]'
    else:
        return r'$p_T^{' + channel_label + '}$ [GeV]'


def extract_hist_name(matrix_name, is_acceptance=True):
    matches = re.findall(r'(\[[^]]*])', matrix_name)
    # case 1: 2 matches, 1D unfold
    if len(matches) == 2:
        remove_string = matches[0]
    # case 2 3 matches, 2D unfold
    elif len(matches) == 4:
        remove_string = matches[2]
        matrix_name = re.sub(r'matrix', 'hist', matrix_name)
    else:
        return
    hist_name = matrix_name.split(remove_string)[0][:-1] + matrix_name.split(remove_string)[1]
    if is_acceptance:
        hist_name = re.sub(r'(__)', '_acceptance__', hist_name)
    else:
        hist_name = re.sub(r'(__)', '_efficiency__', hist_name)
    return hist_name


class ISRAnalyzer(Analyzer):
    def __init__(self, file_pather,
                 signal_dict, bg_dict,
                 channel, period,
                 pt_hist_name, mass_hist_name,
                 pt_matrix_name, mass_matrix_name,
                 mass_edges=None, pt_edges=None,
                 file_path_postfix='',
                 hist_path_postfix=''):

        super().__init__(file_pather, signal_dict, bg_dict)
        self.analysis_name = "ISR"

        self.file_path_postfix = file_path_postfix
        self.hist_path_postfix = hist_path_postfix
        self.hist_producer.set_base_configs(period, channel, file_path_postfix, hist_path_postfix)
        self.plotter.set_period_channel(period, channel, hist_path_postfix)
        self.channel = channel
        self.period = period

        self.pt_hist_name = pt_hist_name
        self.mass_hist_name = mass_hist_name
        self.pt_matrix_name = pt_matrix_name
        self.mass_matrix_name = mass_matrix_name

        self.mass_edges = mass_edges
        self.pt_edges = pt_edges
        self.postfix_mass_window = None
        self.postfix_pt_window = None
        self.is_2d = True
        self.set_postfix_for_1d_case()

        self.mass_unfolders = []
        self.pt_unfolders = []

        # ISR results
        # self.data_isr_hists = None
        # self.data_bg_subtracted__isr_hists
        self.folded_isr_hists = None  # DY fake subtracted
        self.unfolded_isr_hists = None
        self.full_phase_isr_hists = None

        self.set_unfolders()
        self.unfold()
        self.acceptance_corrections()

    def set_postfix_for_1d_case(self):
        if self.mass_edges and self.pt_edges:
            self.is_2d = False

            self.postfix_mass_window = [
                '_dimass_{}to{}'.format(float(self.mass_edges[index][0]),
                                        float(self.mass_edges[index][1]))
                for index in range(len(self.mass_edges))
            ]
            self.postfix_pt_window = [
                '_dipt_{}to{}'.format(float(self.pt_edges[index][0]),
                                      float(self.pt_edges[index][1]))
                for index in range(len(self.pt_edges))
            ]

    def get_pt_mass_hist(self, bin_width_norm=True, is_mass=True, sample_type='data', subtract_bg=False):

        if is_mass:
            hist_name = self.mass_hist_name
            second_axis_postfix = self.postfix_pt_window
        else:
            hist_name = self.pt_hist_name
            second_axis_postfix = self.postfix_mass_window
        self.hist_producer.set_configs(hist_name,
                                       bin_width_norm=bin_width_norm,
                                       second_axis_postfix=second_axis_postfix,
                                       set_text=True)
        if self.is_2d:
            force_1d_output = True
        else:
            force_1d_output = False

        if sample_type == 'data':
            return self.hist_producer.get_data_hist(force_1d_output=force_1d_output,
                                                    bg_subtraction=subtract_bg)
        else:
            return self.hist_producer.get_total_expectation_hist(force_1d_output=force_1d_output,
                                                                 mc_type=sample_type)

    def get_unfold_input_hists_and_matrix(self, unfolder_name="mass"):
        if unfolder_name == "mass":
            hist_name = self.mass_hist_name
            matrix_name = self.mass_matrix_name
            second_axis_edges = self.postfix_pt_window
        else:
            hist_name = self.pt_hist_name
            matrix_name = self.pt_matrix_name
            second_axis_edges = self.postfix_mass_window

        # set histogram config
        self.hist_producer.set_configs(hist_name,
                                       second_axis_postfix=second_axis_edges)
        # set data histogram
        input_hist = self.hist_producer.get_data_hist()  # list of 1D or one 2D
        # set background histograms
        bg_hist = self.hist_producer.get_total_expectation_hist(mc_type='bg')

        # set matrix
        self.hist_producer.set_configs(matrix_name,
                                       second_axis_postfix=second_axis_edges)
        matrix = self.hist_producer.get_total_expectation_hist(mc_type='signal')
        return input_hist, bg_hist, matrix

    def set_unfolder(self, unfolder_name='mass'):
        input_hist, bg_hist, matrix = self.get_unfold_input_hists_and_matrix(unfolder_name)
        if unfolder_name == "mass":
            unfolder_container = self.mass_unfolders
        else:
            unfolder_container = self.pt_unfolders

        if self.is_2d:
            # UnFolder will extract 1D from 2D input
            unfolder_container.append(UnFolder(matrix, input_hist, bg_hist))
        else:
            for i in range(len(input_hist)):
                if bg_hist is None:
                    unfolder_container.append(UnFolder(matrix[i], input_hist[i]))
                else:
                    unfolder_container.append(UnFolder(matrix[i], input_hist[i], bg_hist[i]))

    def set_unfolders(self):
        self.set_unfolder()
        self.set_unfolder(unfolder_name='pt')

    def get_unfolders_config(self):
        return self.channel, self.period, self.hist_path_postfix

    def unfold(self):
        for i in range(len(self.mass_unfolders)):
            self.mass_unfolders[i].do_unfold()

        folded_pt = []
        unfolded_pt = []
        for i in range(len(self.pt_unfolders)):
            self.pt_unfolders[i].do_unfold()

            folded_pt.append(self.pt_unfolders[i].get_folded_hist())
            unfolded_pt.append(self.pt_unfolders[i].get_unfolded_hist())

        unfolded_mass = self.mass_unfolders[0].get_unfolded_hist()
        self.unfolded_isr_hists = ISRHists(unfolded_mass, *unfolded_pt,
                                           mass_windows=self.mass_edges)
        folded_mass = self.mass_unfolders[0].get_folded_hist()
        self.folded_isr_hists = ISRHists(folded_mass, *folded_pt,
                                         mass_windows=self.mass_edges)

    def acceptance_correction(self, unfolder_name='mass'):
        if unfolder_name == 'mass':
            unfolder = self.mass_unfolders
        else:
            unfolder = self.pt_unfolders

        full_phase_data_hist = []
        for i in range(len(unfolder)):
            matrix_name = unfolder[i].response_matrix.hist_name
            full_phase_hist_name = extract_hist_name(matrix_name)
            unfolded_data_hist = unfolder[i].get_unfolded_hist()
            unfolded_data_hist.set_hist_name(full_phase_hist_name)

            unfolded_mc_hist = unfolder[i].get_expectation_hist(use_matrix=True)
            full_phase_hist = self.hist_producer.get_signal_hist(full_phase_hist_name)
            acceptance = full_phase_hist.divide(unfolded_mc_hist)

            full_phase_data_hist.append(unfolded_data_hist.multiply(acceptance))
        return full_phase_data_hist

    def acceptance_corrections(self):
        full_phase_mass_data_hist = self.acceptance_correction()
        full_phase_pt_data_hist = self.acceptance_correction('pt')

        self.full_phase_isr_hists = ISRHists(full_phase_mass_data_hist[0], *full_phase_pt_data_hist,
                                             mass_windows=self.mass_edges)

    def make_mc_isr_hists(self, unfolded=False, use_matrix=True):
        mass_mc_hist = self.mass_unfolders[0].get_expectation_hist(unfolded=unfolded, use_matrix=use_matrix)
        pt_mc_hists = []
        for i in range(len(self.pt_unfolders)):
            pt_mc_hists.append(self.pt_unfolders[i].get_expectation_hist(unfolded=unfolded, use_matrix=use_matrix))
        out_isr_hists = ISRHists(mass_mc_hist, *pt_mc_hists,
                                 mass_windows=self.mass_edges)
        return out_isr_hists

    def get_mc_gen_isr_hists(self, is_full_phase=True):
        mass_hist_name = extract_hist_name(self.mass_unfolders[0].response_matrix.hist_name, is_full_phase)
        self.hist_producer.set_configs(mass_hist_name)
        mc_mass_hist = self.hist_producer.get_total_expectation_hist(mc_type='signal')

        mc_pt_hists = []
        for i in range(len(self.pt_unfolders)):
            pt_hist_name = extract_hist_name(self.pt_unfolders[i].response_matrix.hist_name, is_full_phase)
            self.hist_producer.set_configs(pt_hist_name)
            mc_pt_hists.append(self.hist_producer.get_total_expectation_hist(mc_type='signal'))
        out_isr_hists = ISRHists(mc_mass_hist, *mc_pt_hists,
                                 mass_windows=self.mass_edges)
        return out_isr_hists

    def get_isr_hists(self, level_name='folded'):
        if level_name == 'folded':
            return self.get_folded_isr_hists()
        elif level_name == 'unfolded':
            return self.get_unfolded_isr_hists()
        elif level_name == 'full_phase':
            return self.get_full_phase_isr_hists()
        else:
            print('check ' + level_name)
            return

    def get_mc_isr_hists(self, level_name='folded'):
        if level_name == 'folded':
            return self.make_mc_isr_hists(unfolded=False)
        elif level_name == 'unfolded':
            return self.make_mc_isr_hists(unfolded=True)
        elif level_name == 'efficiency':
            return self.get_mc_gen_isr_hists(False)
        elif level_name == 'full_phase':
            return self.get_mc_gen_isr_hists(True)
        else:
            print('check ' + level_name)
            return

    def get_folded_isr_hists(self):
        return self.folded_isr_hists

    def get_unfolded_isr_hists(self):
        return self.unfolded_isr_hists

    def get_full_phase_isr_hists(self):
        return self.full_phase_isr_hists

    def draw_detector_pt_mass_plots(self, bin_width_norm=True,
                                    is_mass=True, subtract_bg=False):

        data_hists = self.get_pt_mass_hist(bin_width_norm,
                                           is_mass, 'data', subtract_bg)
        mc_type = 'total'
        if subtract_bg:
            mc_type = 'signal'
        total_expectation_hists = self.get_pt_mass_hist(bin_width_norm, is_mass, mc_type)

        for index in range(len(data_hists)):
            self.plotter.add_hist(data_hists[index], color='black', xerr=True, histtype='errorbar')
            as_stack = True
            if subtract_bg:
                as_stack = False
            self.plotter.add_hist(total_expectation_hists[index], is_denominator=True, as_stack=as_stack)

            x_axis_label = make_axis_label(self.channel, is_mass)
            self.draw_data_expectation(x_axis_label)

    def draw_unfolded_plots(self, bin_width_norm=True, draw_mass=True):

        mc_isr_hist = self.get_mc_isr_hists('unfolded')
        data_isr_hist = self.get_isr_hists('unfolded')

        if draw_mass:
            unfolded_data = data_isr_hist.get_mass_hist(bin_width_norm)
            expectation = mc_isr_hist.get_mass_hist(bin_width_norm)

            self.plotter.add_hist(unfolded_data, color='black', xerr=True,
                                  histtype='errorbar')
            self.plotter.add_hist(expectation, is_denominator=True)
            x_axis_label = make_axis_label(self.channel, draw_mass)
            self.draw_data_expectation(x_axis_label)
        else:
            unfolded_data = data_isr_hist.get_pt_hist(-1, bin_width_norm)
            expectation = mc_isr_hist.get_pt_hist(-1, bin_width_norm)

            x_axis_label = make_axis_label(self.channel, draw_mass)
            for index in range(len(unfolded_data)):
                self.plotter.add_hist(unfolded_data[index], color='black', xerr=True,
                                      histtype='errorbar')
                self.plotter.add_hist(expectation[index], is_denominator=True)
                self.draw_data_expectation(x_axis_label)

    def draw_acceptance_or_efficiency_plots(self, draw_mass=True, level='full_phase'):
        unfolded_isr_hist = self.get_mc_isr_hists('unfolded')
        full_phase_isr_hist = self.get_mc_isr_hists(level)

        if level == "full_phase":
            y_axis_label = r"Acc.$\times$ Eff."
            out_name = 'acceptance_efficiency'
        else:
            y_axis_label = "Eff."
            out_name = 'efficiency'

        if draw_mass:
            unfolded = unfolded_isr_hist.get_mass_hist()
            full_phase = full_phase_isr_hist.get_mass_hist()

            hist = unfolded.divide(full_phase)
            x_axis_label = make_axis_label(self.channel)

            self.plotter.add_hist(hist, color='black', xerr=True, histtype='errorbar')
            self.draw(y_axis_label, x_axis_label, out_name)
        else:
            unfolded = unfolded_isr_hist.get_pt_hist(-1)
            full_phase = full_phase_isr_hist.get_pt_hist(-1)

            x_axis_label = make_axis_label(self.channel)
            for index in range(len(unfolded)):
                hist = unfolded[index].divide(full_phase[index])

                self.plotter.add_hist(hist, color='black', xerr=True, histtype='errorbar')
                self.draw(y_axis_label, x_axis_label, out_name)

    def draw_isr_plots(self, level_name, binned_mean=True, show_mc=False):

        # get isr_hists
        isr_hists = self.get_isr_hists(level_name)
        data_mass_df, data_pt_df = isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        label = level_name + ' Data (' + self.channel + ')'
        self.plotter.draw_isr_data_frame(data_mass_df, data_pt_df,
                                         ymin=14, ymax=30,
                                         new_fig=True,
                                         color='black', label=label)
        if show_mc:
            isr_hists = self.get_mc_isr_hists(level_name)
            data_mass_df, data_pt_df = isr_hists.get_isr_dataframe(binned_mean=binned_mean)
            label = level_name + ' MC'
            self.plotter.draw_isr_data_frame(data_mass_df, data_pt_df,
                                             ymin=14, ymax=30,
                                             new_fig=False,
                                             color='gray', label=label)

    def draw_uncertainty_plot(self, level="full_phase"):
        isr_hists = self.get_isr_hists(level)
        mass_df, pt_df = isr_hists.get_isr_dataframe()
        self.plotter.draw_isr_error_df(pt_df)

    def draw_data_expectation(self, x_axis_label=''):
        self.plotter.create_subplots(2, 1, figsize=(10, 10), height_ratios=[1, 0.3])
        self.plotter.set_current_axis(0, 0)
        self.plotter.draw()
        self.plotter.set_y_axis_config(log_scale=True, set_grid=True)
        self.plotter.set_x_axis_config(remove_tick_labels=True, set_grid=True)
        self.plotter.set_axis_config()
        self.plotter.set_experiment_label()

        self.plotter.set_current_axis(1, 0)
        self.plotter.draw_ratio()
        self.plotter.set_y_axis_config(set_min=0.5, set_max=1.5)
        self.plotter.get_axis(1, 0).set_xlabel(x_axis_label, fontsize=30)
        self.plotter.save_fig()

    def draw(self, y_axis_label='', x_axis_label='', out_name=''):
        self.plotter.create_subplots(1, 1, figsize=(10, 10))
        self.plotter.set_current_axis(0, 0)
        self.plotter.draw()
        self.plotter.set_y_axis_config(label=y_axis_label, set_min=0.0, set_max=1.05, set_grid=True)
        self.plotter.set_x_axis_config(label=x_axis_label, set_grid=True)
        # self.plotter.set_axis_config()
        self.plotter.set_experiment_label('Simulation')
        self.plotter.save_fig(out_name)
