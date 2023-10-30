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


class Thesis:
    def __init__(self):
        self.analyzers = {}  # analysis

    def create_isr_analyzer(self, experiment_name, channel, period,
                            pt_hist_name, mass_hist_name,
                            pt_matrix_name, mass_matrix_name,
                            mass_edges=None, pt_edges=None,
                            file_path_postfix='',
                            hist_path_postfix=''):

        if experiment_name == "CMS":
            cms = FilePather("/Users/jhkim/cms_snu/data/Ultralegacy/new_default/",
                             "data_config.json",
                             "mc_config.json",
                             sys_on=False)
            signal_groups = {"Drell-Yan": ["DYJets"]}
            bg_groups = {r"$\tau\tau$:tau": ["DYJets"],
                         "VV": ["ZZ", "WZ", "WW"],
                         r"${t}\bar{t}$": ["TTLL"],
                         "top": ["top"]}
        else:
            return

        analyzer = ISRAnalyzer(cms, signal_groups, bg_groups,
                               channel, period, pt_hist_name, mass_hist_name,
                               pt_matrix_name, mass_matrix_name,
                               mass_edges, pt_edges,
                               file_path_postfix,
                               hist_path_postfix)

        self.analyzers[experiment_name + "_" + channel + period +
                       "_" + file_path_postfix + hist_path_postfix] = analyzer

    def draw_isr_plots(self, *hist_info, ymin=14, ymax=30):

        HistInfo = namedtuple('HistInfo',
                              ['name', 'hist_type', 'level_name', 'binned_mean', 'label', 'kwargs'])
        # FIXME properly set plot path
        plotter = Plotter('CMS', './plots/CMS')

        first_plot = True
        for info in hist_info:
            info = HistInfo(info[0], info[1], info[2], info[3], info[4], info[5])
            if info.hist_type == "data":
                isr_hists = self.analyzers[info.name].get_data_isr_hists(info.level_name)
            else:
                isr_hists = self.analyzers[info.name].get_mc_isr_hists(info.level_name)

            mass_df, pt_df = isr_hists.get_isr_dataframe(binned_mean=info.binned_mean)
            label = info.label
            plotter.draw_isr_data_frame(mass_df, pt_df,
                                        ymin=ymin, ymax=ymax,
                                        new_fig=first_plot,
                                        label=label, **info.kwargs)
            first_plot = False
