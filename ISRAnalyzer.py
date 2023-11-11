from UnFolder import UnFolder
from Analyzer import Analyzer
from ISRHists import ISRHists
from collections import namedtuple
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


def extract_hist_name_from_matrix_name(matrix_name, replacement='acceptance'):
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
    hist_name = re.sub(r'(__)', '_' + replacement + '__', hist_name)
    return hist_name


class ISRAnalyzer(Analyzer):
    def __init__(self, file_pather,
                 signal_dict, bg_dict,
                 channel, period,
                 pt_hist_name, mass_hist_name,
                 pt_matrix_name, mass_matrix_name,
                 mass_edges=None, pt_edges=None,
                 file_path_postfix='',
                 hist_path_postfix='',
                 iterative_unfold=False):

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
        self.do_iterative = iterative_unfold

        self.mass_unfolders = []
        self.pt_unfolders = []
        # TODO method to return number of data points

        # results
        # self.data_isr_hists = None
        # self.data_bg_subtracted__isr_hists
        # levels: folded, unfolded, full_phase
        self.folded_isr_hists = None  # DY fake subtracted
        self.unfolded_isr_hists = None
        self.full_phase_isr_hists = None

        self.analysis_finished = False

    def do_analysis(self):
        self.set_unfolders()
        self.unfold()
        self.acceptance_corrections()

        # generator level study

    def number_of_measurements(self):
        if self.full_phase_isr_hists is None:
            print("Run analysis first...")
            return
        else:
            return self.full_phase_isr_hists.number_of_measurements()

    def channel_label(self):
        if self.channel == "ee":
            return r"$ee$"
        else:
            return r"$\mu\mu$"

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

    def get_pt_mass_hist(self, bin_width_norm=True, normalize=False,
                         is_mass=True, sample_type='data', subtract_bg=False, group_name=None):

        if is_mass:
            hist_name = self.mass_hist_name
            second_axis_postfix = self.postfix_pt_window
        else:
            hist_name = self.pt_hist_name
            second_axis_postfix = self.postfix_mass_window
        self.hist_producer.set_configs(hist_name,
                                       bin_width_norm=bin_width_norm,
                                       second_axis_postfix=second_axis_postfix,
                                       normalize=normalize,
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
                                                                 mc_type=sample_type,
                                                                 target_group=group_name)

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
            unfolder_container.append(UnFolder(matrix, input_hist, bg_hist,
                                               iterative=self.do_iterative))
        else:
            for i in range(len(input_hist)):
                if bg_hist is None:
                    unfolder_container.append(UnFolder(matrix[i], input_hist[i], iterative=self.do_iterative))
                else:
                    bg_norm_test = False  # FIXME
                    unfolder_container.append(UnFolder(matrix[i], input_hist[i], bg_hist[i],
                                                       bg_norm=bg_norm_test, iterative=self.do_iterative))

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

            if not self.do_iterative:
                folded_pt.append(self.pt_unfolders[i].get_folded_hist())  # FIXME how to handle?
            unfolded_pt.append(self.pt_unfolders[i].get_unfolded_hist())

        unfolded_mass = self.mass_unfolders[0].get_unfolded_hist()
        self.unfolded_isr_hists = ISRHists(unfolded_mass, *unfolded_pt,
                                           mass_windows=self.mass_edges,
                                           pt_window=self.pt_edges)
        if not self.do_iterative:
            folded_mass = self.mass_unfolders[0].get_folded_hist()
            self.folded_isr_hists = ISRHists(folded_mass, *folded_pt,
                                             mass_windows=self.mass_edges,
                                             pt_window=self.pt_edges)

    def acceptance_correction(self, unfolder_name='mass'):
        if unfolder_name == 'mass':
            unfolder = self.mass_unfolders
        else:
            unfolder = self.pt_unfolders

        full_phase_data_hist = []
        for i in range(len(unfolder)):
            matrix_name = unfolder[i].response_matrix.hist_name
            full_phase_hist_name = extract_hist_name_from_matrix_name(matrix_name, 'acceptance')
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
                                             mass_windows=self.mass_edges,
                                             pt_window=self.pt_edges)

    def make_mc_isr_hists(self, unfolded=False, use_matrix=True, first_bin=0):
        mass_mc_hist = self.mass_unfolders[0].get_expectation_hist(unfolded=unfolded, use_matrix=use_matrix,
                                                                   first_bin=first_bin)
        pt_mc_hists = []
        for i in range(len(self.pt_unfolders)):
            pt_mc_hists.append(self.pt_unfolders[i].get_expectation_hist(unfolded=unfolded, use_matrix=use_matrix,
                                                                         first_bin=first_bin))
        out_isr_hists = ISRHists(mass_mc_hist, *pt_mc_hists,
                                 mass_windows=self.mass_edges,
                                 pt_window=self.pt_edges)
        return out_isr_hists

    def get_mc_gen_isr_hists(self, target_name):
        mass_hist_name = extract_hist_name_from_matrix_name(self.mass_unfolders[0].response_matrix.hist_name,
                                                            target_name)
        self.hist_producer.set_configs(mass_hist_name)
        mc_mass_hist = self.hist_producer.get_total_expectation_hist(mc_type='signal')

        mc_pt_hists = []
        for i in range(len(self.pt_unfolders)):
            pt_hist_name = extract_hist_name_from_matrix_name(self.pt_unfolders[i].response_matrix.hist_name,
                                                              target_name)
            self.hist_producer.set_configs(pt_hist_name)
            mc_pt_hists.append(self.hist_producer.get_total_expectation_hist(mc_type='signal'))
        out_isr_hists = ISRHists(mc_mass_hist, *mc_pt_hists,
                                 mass_windows=self.mass_edges,
                                 pt_window=self.pt_edges)
        return out_isr_hists

    def get_data_isr_hists(self, level_name='folded'):
        if level_name == 'folded':
            return self.get_folded_isr_hists()
        elif level_name == 'unfolded':
            return self.get_unfolded_isr_hists()
        elif level_name == 'full_phase':
            return self.get_full_phase_isr_hists()
        else:
            print('check ' + level_name)
            return

    def get_mc_isr_hists(self, level_name='folded', first_bin=0):
        if level_name == 'folded':
            return self.make_mc_isr_hists(unfolded=False, first_bin=first_bin)
        elif level_name == 'unfolded':
            return self.make_mc_isr_hists(unfolded=True, first_bin=first_bin)
        elif level_name == 'efficiency':
            return self.get_mc_gen_isr_hists('efficiency')
        elif level_name == 'trigger':
            return self.get_mc_gen_isr_hists('trigger')
        elif level_name == 'full_phase':
            return self.get_mc_gen_isr_hists('acceptance')
        else:
            print('check ' + level_name)
            return

    def get_folded_isr_hists(self):
        return self.folded_isr_hists

    def get_unfolded_isr_hists(self):
        return self.unfolded_isr_hists

    def get_full_phase_isr_hists(self):
        return self.full_phase_isr_hists

    def get_isr_df(self, hist_type, level_name, binned_mean=True):

        if hist_type == 'data':
            isr_hists = self.get_data_isr_hists(level_name)
        else:
            isr_hists = self.get_mc_isr_hists(level_name)

        mass_df, pt_df = isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        return mass_df, pt_df

    def draw_detector_pt_mass_plots(self, bin_width_norm=True,
                                    is_mass=True, subtract_bg=False):

        data_hists = self.get_pt_mass_hist(bin_width_norm=bin_width_norm, is_mass=is_mass,
                                           sample_type='data', subtract_bg=subtract_bg)
        mc_type = 'total'
        if subtract_bg:
            mc_type = 'signal'
        total_expectation_hists = self.get_pt_mass_hist(bin_width_norm=bin_width_norm, is_mass=is_mass,
                                                        sample_type=mc_type)
        x_axis_label = make_axis_label(self.channel, is_mass)
        if is_mass:
            x_fig_size = 10
            y_fig_size = 10
            data_marker_size = 17
        else:
            x_fig_size = 10
            y_fig_size = 10
            data_marker_size = 20

        for index in range(len(data_hists)):
            self.plotter.add_hist(data_hists[index], color='black', xerr=True, histtype='errorbar',
                                  markersize=data_marker_size)
            as_stack = True
            if subtract_bg:
                as_stack = False
            self.plotter.add_hist(total_expectation_hists[index], is_denominator=True, as_stack=as_stack,
                                  linewidth=2)
            self.draw_data_expectation(x_fig_size=x_fig_size, y_fig_size=y_fig_size,
                                       x_axis_label=x_axis_label,
                                       set_y_min=1e-1,
                                       set_ratio_min=0.4, set_ratio_max=1.6)

    def draw_background_ratio(self, is_mass=True):

        data_hists = self.get_pt_mass_hist(is_mass=is_mass, sample_type='data')
        total_bg_hists = self.get_pt_mass_hist(is_mass=is_mass, sample_type='bg')

        x_fig_size = 10

        bg_hists = {}
        for bg_group in self.bg_dict:
            bg_group_label = bg_group.split(":")[0]
            bg_hists[bg_group_label] = (self.get_pt_mass_hist(is_mass=is_mass, sample_type='bg', group_name=bg_group))
        x_axis_label = make_axis_label(self.channel, is_mass=is_mass)

        for index in range(len(data_hists)):
            total_ratio = total_bg_hists[index].divide(data_hists[index])

            if is_mass:
                text = r"$\mathit{p_T } < $" + "{:.0f}".format(float(self.pt_edges[index][1])) + " GeV"
            else:
                text = ("{:.0f}".format(float(self.mass_edges[index][0])) + r"$ < \mathit{m} < $" +
                        "{:.0f}".format(float(self.mass_edges[index][1])) + " GeV")

            self.plotter.set_text(text=text)
            self.plotter.add_hist(total_ratio, xerr=True, histtype='step', label="Total bg", linewidth=2)

            for label, bg_hist in bg_hists.items():
                bg_ratio = bg_hist[index].divide(data_hists[index])
                self.plotter.add_hist(bg_ratio, xerr=True, histtype='errorbar', label=label,
                                      markersize=20)

            out_name = "bg_fraction"
            if is_mass:
                out_name += "_mass"
            else:
                out_name += "_pt"

            self.draw(x_fig_size=x_fig_size, y_axis_label="Fraction", legend_loc='lower right',
                      do_mpl_magic=False,
                      x_axis_label=x_axis_label, set_min=1e-5, set_max=1, set_y_log=True,
                      out_name=out_name)

    def draw_simulation_pt_mass_plots(self, level, draw_mass=True, bin_width_norm=True):
        isr_hists = self.get_mc_isr_hists(level)
        x_axis_label = make_axis_label(self.channel, draw_mass)
        if draw_mass:
            hist = isr_hists.get_mass_hist(bin_width_norm)
            self.plotter.add_hist(hist, color='black', xerr=True, histtype='step', linewidth=2)
            self.draw(x_fig_size=10, x_axis_label=x_axis_label, set_y_log=True)
        else:
            index_to_draw = [0, 2, 4, 6, 7]
            hists = isr_hists.get_pt_hist(-1, bin_width_norm=True, normalize=True)
            self.plotter.use_text_as_label = True
            for index in range(len(hists)):
                if index not in index_to_draw:
                    continue
                is_denominator = False
                if index == 4:
                    is_denominator = True
                self.plotter.add_hist(hists[index], is_denominator=is_denominator,
                                      xerr=True, histtype='step', linewidth=2)
                # self.draw(x_axis_label=x_axis_label)
            self.draw_data_expectation(plot_label="Simulation", set_y_min=1e-4, set_y_max=0.5,
                                       x_axis_label=x_axis_label, set_ratio_min=0.0, set_ratio_max=2,
                                       do_mpl_magic=False, legend_loc='lower left',
                                       ratio_label='Ratio', set_ratio_log=False)
            # self.draw(x_axis_label=x_axis_label, set_y_log=True, set_legend=False, do_mpl_magic=False,
            #           set_x_log=False)

    def draw_unfolded_plots(self, bin_width_norm=True, draw_mass=True):

        mc_isr_hist = self.get_mc_isr_hists('unfolded')
        data_isr_hist = self.get_data_isr_hists('unfolded')

        x_axis_label = make_axis_label(self.channel, draw_mass)
        if draw_mass:
            unfolded_data = data_isr_hist.get_mass_hist(bin_width_norm)
            expectation = mc_isr_hist.get_mass_hist(bin_width_norm)

            self.plotter.add_hist(unfolded_data, color='black', xerr=True,
                                  histtype='errorbar')
            self.plotter.add_hist(expectation, is_denominator=True)
            self.draw_data_expectation(x_axis_label=x_axis_label)
        else:
            unfolded_data = data_isr_hist.get_pt_hist(-1, bin_width_norm)
            expectation = mc_isr_hist.get_pt_hist(-1, bin_width_norm)

            for index in range(len(unfolded_data)):
                self.plotter.add_hist(unfolded_data[index], color='black', xerr=True,
                                      histtype='errorbar')
                self.plotter.add_hist(expectation[index], is_denominator=True)
                self.draw_data_expectation(x_axis_label=x_axis_label)

    def draw_acceptance_or_efficiency_plots(self, *hist_types, draw_mass=True):
        # event acc x eff: reco & gen / gen
        # event eff: reco & gen / gen with acceptance cut
        # SF: reco & gen with SF applied/ reco & gen without SF
        first_hist = True
        hists = []
        out_name = ''
        for hist_type in hist_types:
            if hist_type == "acc_eff":
                nominator_isr_hist = self.get_mc_isr_hists('unfolded')
                denominator_isr_hist = self.get_mc_isr_hists('full_phase')
                y_axis_label = r"Acc.$\times$ Eff."
                out_name += 'acceptance_efficiency'
            elif hist_type == "acc":
                nominator_isr_hist = self.get_mc_isr_hists('efficiency')
                denominator_isr_hist = self.get_mc_isr_hists('full_phase')
                y_axis_label = r"Acceptance"
                out_name += 'acceptance'
            elif hist_type == "eff":
                nominator_isr_hist = self.get_mc_isr_hists('unfolded')
                denominator_isr_hist = self.get_mc_isr_hists('efficiency')
                y_axis_label = "Eff."
                out_name += 'efficiency'
            elif hist_type == "trig":
                nominator_isr_hist = self.get_mc_isr_hists('trigger')
                denominator_isr_hist = self.get_mc_isr_hists('efficiency')
                y_axis_label = "Trigger eff."
                out_name += 'trigger'
            else:
                nominator_isr_hist = self.get_mc_isr_hists('unfolded', first_bin=1)
                denominator_isr_hist = self.get_mc_isr_hists('unfolded')
                y_axis_label = "SF"
                out_name += 'SF'

            if draw_mass:
                nominator = nominator_isr_hist.get_mass_hist()
                denominator = denominator_isr_hist.get_mass_hist()

                hist = nominator.divide(denominator)
                if first_hist:
                    hists.append({})
                    first_hist = False
                hists[0].update({y_axis_label: hist})
            else:
                nominator = nominator_isr_hist.get_pt_hist(-1)
                denominator = denominator_isr_hist.get_pt_hist(-1)

                for index in range(len(nominator)):
                    if first_hist:
                        hists.append({})
                    hist = nominator[index].divide(denominator[index])
                    hists[index].update({y_axis_label: hist})
                first_hist = False

        x_axis_label = make_axis_label(self.channel, is_mass=draw_mass)
        if draw_mass:
            for label, hist in hists[0].items():
                self.plotter.add_hist(hist,
                                      xerr=True, histtype='errorbar', label=label, markersize=20)
            self.draw(y_axis_label="Fraction", x_axis_label=x_axis_label,
                      set_min=0.0, set_max=1.05, out_name=out_name)
        else:
            for i in range(len(hists)):
                for label, hist in hists[i].items():
                    self.plotter.add_hist(hist, xerr=True, histtype='errorbar', label=label, markersize=20)
                self.draw(y_axis_label="Fraction", x_axis_label=x_axis_label,
                          set_min=0.0, set_max=1.05, out_name=out_name)

    def draw_isr_plots(self, *hist_info, ymin=14, ymax=30):
        # (hist_type, leve_name, binned_mean, label, kwargs)
        # get isr_hists
        HistInfo = namedtuple('HistInfo',
                              ['hist_type', 'level_name', 'binned_mean', 'label', 'kwargs'])
        first_plot = True
        for info in hist_info:
            info = HistInfo(info[0], info[1], info[2], info[3], info[4])
            if info.hist_type == "data":
                isr_hists = self.get_data_isr_hists(info.level_name)
                label = "Data " + self.channel_label()
            else:
                isr_hists = self.get_mc_isr_hists(info.level_name)
                label = ""

            mass_df, pt_df = isr_hists.get_isr_dataframe(binned_mean=info.binned_mean)
            label = label + info.label
            self.plotter.draw_isr_data_frame(mass_df, pt_df,
                                             ymin=ymin, ymax=ymax,
                                             new_fig=first_plot,
                                             label=label, **info.kwargs)
            first_plot = False

    def draw_uncertainty_plot(self, level="full_phase"):
        isr_hists = self.get_data_isr_hists(level)
        mass_df, pt_df = isr_hists.get_isr_dataframe()
        self.plotter.draw_isr_error_df(pt_df)

    def draw_data_expectation(self, x_fig_size=10, y_fig_size=10, x_axis_label='',
                              plot_label="Preliminary", do_mpl_magic=True, legend_loc='best',
                              ratio_label=None, set_ratio_log=False,
                              set_y_min=None, set_y_max=None,
                              set_ratio_min=None, set_ratio_max=None):

        self.plotter.create_subplots(2, 1, figsize=(x_fig_size, y_fig_size), height_ratios=[1, 0.3])
        self.plotter.set_current_axis(0, 0)
        self.plotter.draw()
        self.plotter.set_y_axis_config(log_scale=True, set_min=set_y_min, set_max=set_y_max)
        self.plotter.set_x_axis_config(remove_tick_labels=True)
        self.plotter.set_axis_config(set_legend=True, do_mpl_magic=do_mpl_magic, legend_loc=legend_loc)
        self.plotter.set_experiment_label(plot_label)

        self.plotter.set_current_axis(1, 0)
        self.plotter.draw_ratio(y_axis_label=ratio_label)
        self.plotter.set_y_axis_config(set_min=set_ratio_min, set_max=set_ratio_max,
                                       set_grid=True, log_scale=set_ratio_log)
        self.plotter.set_x_axis_config(label=x_axis_label)
        self.plotter.save_fig()

    def draw(self, x_fig_size=10, y_axis_label='', x_axis_label='', do_mpl_magic=True, set_legend=True,
             legend_loc='best', set_y_log=False, set_x_log=False,
             set_min=None, set_max=None, out_name=''):
        self.plotter.create_subplots(1, 1, figsize=(x_fig_size, 10))
        self.plotter.set_current_axis(0, 0)
        self.plotter.draw()
        self.plotter.set_y_axis_config(label=y_axis_label, log_scale=set_y_log,
                                       set_min=set_min, set_max=set_max, set_grid=True)
        self.plotter.set_x_axis_config(label=x_axis_label, set_grid=True, log_scale=set_x_log)
        self.plotter.set_axis_config(legend_loc=legend_loc, do_mpl_magic=do_mpl_magic, set_legend=set_legend)
        self.plotter.set_experiment_label('Simulation')
        self.plotter.save_fig(out_name)
