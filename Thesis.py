from FilePather import FilePather
from ISRAnalyzer import ISRAnalyzer
from collections import namedtuple
from Plotter import Plotter


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


def draw_comparison_plot(plotter, y_axis_label='', x_axis_label='', legend_loc='best',
                         ratio_label=None, ratio_y_min=0.5, ratio_y_max=1.5,
                         out_name=''):

    plotter.create_subplots(2, 1, figsize=(10, 10), height_ratios=[1, 0.3])
    plotter.set_current_axis(0, 0)
    plotter.draw()
    plotter.set_y_axis_config(label=y_axis_label, log_scale=True, set_grid=True)
    plotter.set_x_axis_config(remove_tick_labels=True, set_grid=True)
    plotter.set_axis_config(legend_loc=legend_loc)
    plotter.set_experiment_label()

    plotter.set_current_axis(1, 0)
    plotter.draw_ratio(y_axis_label=ratio_label)
    plotter.set_y_axis_config(set_min=ratio_y_min, set_max=ratio_y_max, set_grid=True)
    plotter.set_x_axis_config(label=x_axis_label)
    plotter.save_fig(out_name=out_name)


class Thesis:
    def __init__(self):
        self.analyzers = {}  # analysis

    def create_isr_analyzer(self, experiment_name, channel, period,
                            pt_hist_name, mass_hist_name,
                            pt_matrix_name, mass_matrix_name,
                            mass_edges=None, pt_edges=None,
                            file_path_postfix='',
                            hist_path_postfix='',
                            use_minnlo=True,
                            do_iterative=False):

        if experiment_name == "CMS":
            cms = FilePather("/Users/jhkim/cms_snu/data/Ultralegacy/default/",
                             "data_config.json",
                             "mc_config.json",
                             sys_on=False)
            signal_prefix = ""
            dy_tau_sample_name = "DYJets"
            if use_minnlo:
                dy_tau_sample_name = "DYJetsToTauTau_MiNNLO"
                if channel == "ee":
                    signal_groups = {"DY (MiNNLO)": ["DYJetsToEE_MiNNLO"]}
                else:
                    signal_groups = {"DY (MiNNLO)": ["DYJetsToMuMu_MiNNLO"]}
            else:
                signal_prefix = "_aMCNLO"
                signal_groups = {"Drell-Yan": ["DYJets"]}

            bg_groups = {r"DY ($\tau\tau$):tau": [dy_tau_sample_name],
                         r"$VV$": ["ZZ", "WZ", "WW"],
                         r"${t}\bar{t}+tW+\bar{t}W$": ["TTLL", "top", "antitop"]}
        else:
            return

        analyzer = ISRAnalyzer(cms, signal_groups, bg_groups,
                               channel, period, pt_hist_name, mass_hist_name,
                               pt_matrix_name, mass_matrix_name,
                               mass_edges, pt_edges,
                               file_path_postfix,
                               hist_path_postfix,
                               iterative_unfold=do_iterative)

        postfix = ''
        if file_path_postfix != "":
            postfix += "_" + file_path_postfix
        if hist_path_postfix != "":
            postfix += "_" + hist_path_postfix

        self.analyzers[experiment_name + "_" + channel + period +
                       postfix + signal_prefix] = analyzer

    def get_combined_results(self, *hist_info):
        # loop over pairs of pt, mass values
        # loop over hists
        # combine
        pass

    def draw_pt_distributions(self, *hist_info, channel="ee",
                              ratio_label=None,
                              ratio_ymin=0.5, ratio_ymax=1.5):
        # hist_info
        # name, hist_type, bin_width_norm, label, kwargs
        HistInfo = namedtuple('HistInfo',
                              ['name', 'hist_type', 'subtract_bg', 'normalize',
                               'label', 'bin_width', 'kwargs'])
        x_axis_label = make_axis_label(channel, False)
        first_hist_info = True
        plotters = []
        for info in hist_info:
            info = HistInfo(info[0], info[1], info[2], info[3], info[4], info[5], info[6])
            hists = self.analyzers[info.name].get_pt_mass_hist(bin_width_norm=info.bin_width,
                                                               is_mass=False,
                                                               sample_type=info.hist_type,
                                                               subtract_bg=info.subtract_bg,
                                                               normalize=info.normalize)
            for index in range(len(hists)):
                if first_hist_info:
                    plotters.append(Plotter('CMS', './plots/CMS'))
                    plotters[index].add_hist(hists[index],
                                             is_denominator=True, **info.kwargs)
                else:
                    plotters[index].add_hist(hists[index],
                                             is_denominator=False, **info.kwargs)
            first_hist_info = False
        for index, plotter in enumerate(plotters):
            draw_comparison_plot(plotter, x_axis_label=x_axis_label,
                                 ratio_label=ratio_label,
                                 ratio_y_min=ratio_ymin, ratio_y_max=ratio_ymax)

    def draw_background_ratio(self, *hist_info, channel="ee",
                              ratio_label=None,
                              ratio_ymin=0.5, ratio_ymax=1.5):

        HistInfo = namedtuple('HistInfo',
                              ['name', 'target_bg', 'label', 'kwargs'])
        x_axis_label = make_axis_label(channel, False)
        first_hist_info = True
        plotters = []

        y_axis_label='Fraction'
        out_name = "bg_ratio"
        for info in hist_info:
            info = HistInfo(info[0], info[1], info[2], info[3])
            data_hists = self.analyzers[info.name].get_pt_mass_hist(is_mass=False, sample_type='data')
            hists = self.analyzers[info.name].get_pt_mass_hist(is_mass=False,
                                                               sample_type='bg', group_name=info[1])

            for index in range(len(data_hists)):
                fraction = hists[index].divide(data_hists[index])
                if first_hist_info:
                    text = ("{:.0f}".format(float(self.analyzers[info.name].mass_edges[index][0])) + r"$ < \mathit{m} < $" +
                            "{:.0f}".format(float(self.analyzers[info.name].mass_edges[index][1])) + " GeV")

                    plotters.append(Plotter('CMS', './plots/CMS'))
                    plotters[index].set_text(text=text)
                    plotters[index].add_hist(fraction, is_denominator=True, **info.kwargs)
                else:
                    plotters[index].add_hist(fraction, is_denominator=False, **info.kwargs)
            first_hist_info = False
        for index, plotter in enumerate(plotters):
            draw_comparison_plot(plotter, y_axis_label, x_axis_label, legend_loc='upper right',
                                 ratio_label=ratio_label,
                                 ratio_y_min=ratio_ymin, ratio_y_max=ratio_ymax,
                                 out_name=out_name)

    def draw_isr_plots(self, *hist_info, ymin=14, ymax=30):

        HistInfo = namedtuple('HistInfo',
                              ['name', 'hist_type', 'level_name', 'binned_mean', 'label', 'kwargs'])
        # FIXME properly set plot path
        plotter = Plotter('CMS', './plots/CMS')

        first_plot = True
        for info in hist_info:
            info = HistInfo(info[0], info[1], info[2], info[3], info[4], info[5])
            mass_df, pt_df = self.analyzers[info.name].get_isr_df(info.hist_type, info.level_name,
                                                                  info.binned_mean)
            label = info.label
            plotter.draw_isr_data_frame(mass_df, pt_df,
                                        ymin=ymin, ymax=ymax,
                                        new_fig=first_plot,
                                        label=label, **info.kwargs)
            first_plot = False
