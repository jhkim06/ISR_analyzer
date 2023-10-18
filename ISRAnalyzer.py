from Analyzer import Analyzer
from Unfolder import UnFolder
from ISRHists import ISRHists
from HistProducer import HistProducer
import re


def get_full_phase_hist_name(hist):
    hist_name = re.sub(r'(__)', '_acceptance__', hist.hist_name)
    return hist_name


class ISRAnalyzer(Analyzer):
    def __init__(self, file_pather, signal_dict=None, bg_dict=None):

        super().__init__(file_pather, signal_dict, bg_dict)
        self.analysis_name = "ISR"

        # major histogram name
        self.pt_hist_name = ''
        self.mass_hist_name = ''
        self.pt_matrix_name = ''
        self.mass_matrix_name = ''

        self.mass_edges = None
        self.pt_edges = None
        self.mass_postfix_for_pt_hist = None
        self.pt_postfix_for_mass_hist = None
        self.use_2D = True  # if 2D, dipt histograms for each mass window in one 1D histogram

        self.input_data_isr_hists = None  # before bkg. subtraction
        self.folded_isr_hists = None
        self.unfolded_isr_hists = None
        self.full_phase_isr_hists = None

        self.mass_unfolder = []
        self.pt_unfolder = []

        self.bg_used = True
        if bg_dict is None:
            self.bg_used = False

    def set_mass_pt_hist_names(self, mass_hist_name, pt_hist_name, mass_matrix_name, pt_matrix_name,
                               mass_edges=None, pt_edges=None):
        self.mass_hist_name = mass_hist_name
        self.mass_matrix_name = mass_matrix_name
        self.pt_hist_name = pt_hist_name
        self.pt_matrix_name = pt_matrix_name

        if mass_edges is not None and pt_edges is not None:
            self.use_2D = False
            self.mass_edges = mass_edges
            self.pt_edges = pt_edges

            self.mass_postfix_for_pt_hist = ['_dimass_{}to{}'.format(
                float(mass_edges[index][0]), float(mass_edges[index][1])) for index in range(len(mass_edges))]
            self.pt_postfix_for_mass_hist = ['_dipt_{}to{}'.format(
                float(pt_edges[index][0]), float(pt_edges[index][1])) for index in range(len(pt_edges))]
        else:
            self.use_2D = True

    def draw_plot(self, x_axis_label):
        # for 2D
        if self.use_2D:
            data_hists = self.hist_producer.get_data_hist().get_1d_hists()
            total_expectation_hists = self.hist_producer.get_total_expectation_hist().get_1d_hists()
        else:
            data_hists = self.hist_producer.get_data_hist()
            total_expectation_hists = self.hist_producer.get_total_expectation_hist()

        for index in range(len(data_hists)):
            self.plotter.draw_comparison(total_expectation_hists[index],
                                         data_hists[index],
                                         x_axis_label=x_axis_label, logy=True)

    # draw pt and mass plots
    def draw_detector_plots(self):
        self.hist_producer.set_configs(self.pt_hist_name,
                                       bin_width_norm=True,
                                       second_axis_postfix=self.mass_postfix_for_pt_hist,
                                       set_text=True)
        x_axis_label = '$\mathit{p_T ^{' + self.get_channel_label() + '}}$ [GeV]'
        self.draw_plot(x_axis_label=x_axis_label)

        # mass
        self.hist_producer.set_configs(self.mass_hist_name,
                                       bin_width_norm=True,
                                       second_axis_postfix=self.pt_postfix_for_mass_hist,
                                       set_text=True)
        x_axis_label = '$\mathit{m ^{' + self.get_channel_label() + '}}$ [GeV]'
        self.draw_plot(x_axis_label=x_axis_label)

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
        if self.bg_used:
            bg_hist = self.hist_producer.get_total_expectation_hist(exp_type='bg')
        else:
            bg_hist = None
        # set matrix
        self.hist_producer.set_configs(matrix_name,
                                       second_axis_postfix=second_axis_edges)
        matrix = self.hist_producer.get_total_expectation_hist(exp_type='signal')
        # TODO return UnFolder()
        return input_hist, bg_hist, matrix

    def set_unfold(self):
        # initialize the list of mass and pt unfolder before append
        self.pt_unfolder.clear()
        self.mass_unfolder.clear()

        mass_data_hist, mass_bg_hist, mass_matrix = self.get_unfold_input_hists_and_matrix()
        pt_data_hist, pt_bg_hist, pt_matrix = self.get_unfold_input_hists_and_matrix(is_mass=False)
        # use default unfold setting
        if self.use_2D:
            self.mass_unfolder.append(UnFolder(mass_matrix, mass_data_hist, mass_bg_hist))
            self.pt_unfolder.append(UnFolder(pt_matrix, pt_data_hist, pt_bg_hist))
            self.input_data_isr_hists = ISRHists(mass_data_hist, pt_data_hist, mass_windows=self.mass_edges)
        else:
            for i in range(len(mass_data_hist)):
                if self.bg_used:
                    self.mass_unfolder.append(UnFolder(mass_matrix[i], mass_data_hist[i], mass_bg_hist[i]))
                else:
                    self.mass_unfolder.append(UnFolder(mass_matrix[i], mass_data_hist[i]))
            for i in range(len(pt_data_hist)):
                if self.bg_used:
                    self.pt_unfolder.append(UnFolder(pt_matrix[i], pt_data_hist[i], pt_bg_hist[i]))
                else:
                    self.pt_unfolder.append(UnFolder(pt_matrix[i], pt_data_hist[i]))
            self.input_data_isr_hists = ISRHists(mass_data_hist[0], *pt_data_hist, mass_windows=self.mass_edges)

    def unfold(self):
        for i in range(len(self.mass_unfolder)):
            self.mass_unfolder[i].do_unfold()

        unfolded_pt = []
        folded_pt = []
        for i in range(len(self.pt_unfolder)):
            self.pt_unfolder[i].do_unfold()

            unfolded_pt.append(self.pt_unfolder[i].get_unfolded_hist())
            folded_pt.append(self.pt_unfolder[i].get_folded_hist())

        unfolded_mass = self.mass_unfolder[0].get_unfolded_hist()
        folded_mass = self.mass_unfolder[0].get_folded_hist()
        self.folded_isr_hists = ISRHists(folded_mass, *folded_pt, mass_windows=self.mass_edges)
        self.unfolded_isr_hists = ISRHists(unfolded_mass, *unfolded_pt, mass_windows=self.mass_edges)

    def acceptance_correction(self, unfolder):
        full_phase_data_hist = []
        for i in range(len(unfolder)):
            unfolded_data_hist = unfolder[i].get_unfolded_hist()
            full_phase_hist_name = get_full_phase_hist_name(unfolded_data_hist)
            unfolded_data_hist.set_hist_name(full_phase_hist_name)

            unfolded_mc_hist = unfolder[i].get_expectation_hist(use_matrix=True)

            local_producer = HistProducer(self.file_pather, self.signal_dict, self.bg_dict)
            local_producer.set_base_configs(self.period, self.channel, self.file_path_postfix)
            full_phase_hist = local_producer.get_signal_hist(full_phase_hist_name)

            acceptance = full_phase_hist.divide(unfolded_mc_hist)
            full_phase_data_hist.append(unfolded_data_hist.multiply(acceptance))
        return full_phase_data_hist

    def acceptance_corrections(self):
        # get full_phase hist name and
        full_phase_mass_data_hist = self.acceptance_correction(self.mass_unfolder)
        full_phase_pt_data_hist = self.acceptance_correction(self.pt_unfolder)

        self.full_phase_isr_hists = ISRHists(full_phase_mass_data_hist[0], *full_phase_pt_data_hist,
                                             mass_windows=self.mass_edges)

    def get_signal_isr_hist(self, hist_name, unfolded=True, is_full_phase=False, unfolder=None):
        local_producer = HistProducer(self.file_pather, self.signal_dict, self.bg_dict)
        local_producer.set_base_configs(self.period, self.channel, self.file_path_postfix)
        if is_full_phase:
            if self.use_2D:  # when 2D, ignore the attached names
                matches = re.findall(r'(\[[^]]*])', hist_name)
                hist_name = matches[0] + "_" + matches[1] + "_" + matches[2]
            local_producer.set_configs(hist_name)
            mc_hist = local_producer.get_total_expectation_hist(exp_type='signal')
        else:
            # mc_hist = self.mass_unfolder[0].get_expectation_hist(unfolded=unfolded, use_matrix=True)
            mc_hist = unfolder.get_expectation_hist(unfolded=unfolded, use_matrix=True)
        return mc_hist

    def make_expectation_isr_hists(self, isr_hists, unfolded=True, is_full_phase=False):
        # get pt/mass hist from data isr_hist and make expectation isr hist
        mass_hist = isr_hists.get_mass_hist()
        pt_hists = isr_hists.get_pt_hist(-1)

        mass_hist_name = mass_hist.hist_name
        mass_mc_hist = self.get_signal_isr_hist(mass_hist_name, unfolded, is_full_phase, self.mass_unfolder[0])

        pt_mc_hists = []
        for i in range(len(self.pt_unfolder)):
            pt_hist_name = pt_hists[i].hist_name
            pt_mc_hist = self.get_signal_isr_hist(pt_hist_name, unfolded, is_full_phase, self.pt_unfolder[i])
            pt_mc_hists.append(pt_mc_hist)
        out_isr_hists = ISRHists(mass_mc_hist, *pt_mc_hists,
                                 mass_windows=self.mass_edges)
        return out_isr_hists

    # draw_isr_plots(self, ("mc", "folded"), ("da")
    def draw_isr_result(self, level_name='folded', color='black', update=False, show_mc=True,
                        binned_mean=False, ymin=14, ymax=30, xmin=50, xmax=300):

        if level_name == "detector":
            isr_hists = self.input_data_isr_hists
            data_mass_df, data_pt_df = self.input_data_isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        elif level_name == 'folded':
            isr_hists = self.folded_isr_hists
            data_mass_df, data_pt_df = self.folded_isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        elif level_name == 'unfolded':
            isr_hists = self.unfolded_isr_hists
            data_mass_df, data_pt_df = self.unfolded_isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        else:
            isr_hists = self.full_phase_isr_hists
            data_mass_df, data_pt_df = self.full_phase_isr_hists.get_isr_dataframe(binned_mean=binned_mean)
        self.plotter.draw_isr_data_frame(data_mass_df, data_pt_df,
                                         ymin=ymin, ymax=ymax,
                                         new_fig= (not update),
                                         color=color, label=level_name + ' data (' + self.get_channel_label() + ')')
        if show_mc:
            unfolded = True
            is_full_phase = False
            if level_name == 'folded':
                unfolded = False
            if level_name == 'full_phase':
                is_full_phase = True
            mc_isr_hists = self.make_expectation_isr_hists(isr_hists,
                                                           unfolded=unfolded, is_full_phase=is_full_phase)
            if level_name == "detector":
                hist_name = self.mass_hist_name
                second_axis_edges = self.pt_postfix_for_mass_hist
                self.hist_producer.set_configs(hist_name,
                                               second_axis_postfix=second_axis_edges)
                mc_mass_hist = self.hist_producer.get_total_expectation_hist(exp_type='total')

                hist_name = self.pt_hist_name
                second_axis_edges = self.mass_postfix_for_pt_hist
                self.hist_producer.set_configs(hist_name,
                                               second_axis_postfix=second_axis_edges)
                mc_pt_hist = self.hist_producer.get_total_expectation_hist(exp_type='total')
                mc_isr_hists = ISRHists(mc_mass_hist[0], *mc_pt_hist,
                                        mass_windows=self.mass_edges)

            mc_mass_df, mc_pt_df = mc_isr_hists.get_isr_dataframe(binned_mean=binned_mean)
            self.plotter.draw_isr_data_frame(mc_mass_df, mc_pt_df,
                                             ymin=ymin, ymax=ymax,
                                             new_fig=False,
                                             mfc='none',
                                             color=color, label=level_name + ' mc (' + self.get_channel_label() + ')')

    def draw_signal_projected_hist(self):
        if self.use_2D:
            projected_hist = self.pt_unfolder[0].get_projected_hist(first_bin=1).get_1d_hists()
            expectation_hist = self.pt_unfolder[0].get_expectation_hist().get_1d_hists()
        else:
            projected_hist = []
            expectation_hist = []
            for i in range(len(self.pt_unfolder)):
                projected_hist.append(self.pt_unfolder[i].get_projected_hist(first_bin=1))
                expectation_hist.append(self.pt_unfolder[i].get_expectation_hist())

        for index in range(len(projected_hist)):
            self.plotter.draw_comparison(expectation_hist[index], projected_hist[index])

    def draw_unfolded_data(self, bin_width_norm=True):
        # make ISRHists for expectation
        mass_mc_hist = []
        pt_mc_hist = []
        for i in range(len(self.mass_unfolder)):
            mass_mc_hist.append(self.mass_unfolder[i].get_expectation_hist(use_matrix=True))
        for i in range(len(self.pt_unfolder)):
            pt_mc_hist.append(self.pt_unfolder[i].get_expectation_hist(use_matrix=True))
        expectation_isr_hist = ISRHists(mass_mc_hist[0], *pt_mc_hist, mass_windows=self.mass_edges)

        unfolded_data = self.unfolded_isr_hists.get_pt_hist(-1, bin_width_norm)
        expectation = expectation_isr_hist.get_pt_hist(-1, bin_width_norm)

        # pt histogram
        for index in range(len(unfolded_data)):
            self.plotter.draw_comparison(expectation[index], unfolded_data[index],
                                         x_axis_label='$\mathit{p_T ^{' + self.get_channel_label() + '}}$ [GeV]')

        # mass histogram
        unfolded_data = self.unfolded_isr_hists.get_mass_hist(bin_width_norm)
        expectation = expectation_isr_hist.get_mass_hist(bin_width_norm)
        self.plotter.draw_comparison(expectation, unfolded_data,
                                     x_axis_label='$\mathit{m ^{' + self.get_channel_label() + '}}$ [GeV]')
        
    def draw_uncertainty_plot(self):
        mass_df, pt_df = self.unfolded_isr_hists.get_isr_dataframe()
        self.plotter.draw_isr_error_df(pt_df)
