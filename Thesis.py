from FilePather import FilePather
from ISRAnalyzer import ISRAnalyzer


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

    def draw_isr_plot_comparisons(self):
        pass
