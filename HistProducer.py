import re


# create Hist from input ROOT files
class HistProducer:

    def __init__(self, experiment_file_pather, signal_dict=None, bg_dict=None):

        self.experiment_file_pather = experiment_file_pather
        self.signal_dict = signal_dict
        self.bg_dict = bg_dict
        self.expectation_dict = None

        if self.signal_dict is not None and self.bg_dict is not None:
            self.expectation_dict = {**self.signal_dict, **self.bg_dict}
        if self.signal_dict is not None and self.bg_dict is None:
            self.expectation_dict = self.signal_dict

        self.period = ''
        self.channel = ''
        self.hist_name = ''
        self.second_axis_postfix = None  # FIXME naming ["_55.0to64.0", ]

        self.file_path_postfix = ''
        self.hist_path_postfix = ''
        self.bin_width_norm = False
        self.normalize = False
        self.axis_steering = ''
        self.set_text_for_plot = False

    def set_base_configs(self, period, channel,
                         file_path_postfix='',
                         hist_path_postfix=''):

        self.period = period
        self.channel = channel
        self.file_path_postfix = file_path_postfix
        self.hist_path_postfix = hist_path_postfix

    def set_configs(self, hist_name, **kwargs):
        self.hist_name = hist_name

        if 'bin_width_norm' in kwargs:
            self.bin_width_norm = kwargs['bin_width_norm']
        else:
            self.bin_width_norm = False

        if 'normalize' in kwargs:
            self.normalize = kwargs['normalize']
        else:
            self.normalize = False

        if 'axis_steering' in kwargs:
            self.axis_steering = kwargs['axis_steering']
        else:
            self.axis_steering = ''

        if 'second_axis_postfix' in kwargs:
            self.second_axis_postfix = kwargs['second_axis_postfix']
        else:
            self.second_axis_postfix = None

        if 'set_text' in kwargs:
            self.set_text_for_plot = kwargs['set_text']
        else:
            self.set_text_for_plot = False

    def get_signal_hist(self, name):

        signal_label = [*self.signal_dict.keys()][0]
        signal_group = self.experiment_file_pather.make_mc_group(self.period, self.channel,
                                                                 *self.signal_dict[signal_label],
                                                                 group_name=signal_label,
                                                                 file_path_postfix=self.file_path_postfix,
                                                                 hist_path_postfix=self.hist_path_postfix)
        return signal_group.get_hist(name)

    def get_data_hist(self, force_1d_output=False, bg_subtracted=False):
        # allow to return a list of histograms
        data = self.experiment_file_pather.make_data_group(self.period, self.channel,
                                                           file_path_postfix=self.file_path_postfix,
                                                           hist_path_postfix=self.hist_path_postfix)
        if self.second_axis_postfix is None:
            data_hist = data.get_hist(self.hist_name)
            data_hist.set_hist_config(self.bin_width_norm, self.axis_steering, self.normalize)
            if force_1d_output:
                return data_hist.get_1d_hists()
            else:
                return data_hist
        else:
            data_hists = []
            for postfix in self.second_axis_postfix:
                data_hist = data.get_hist(self.hist_name + postfix)
                data_hist.set_hist_config(self.bin_width_norm, self.axis_steering, self.normalize)
                data_hists.append(data_hist)
                if self.set_text_for_plot:
                    edges = re.findall(r'\d+\.\d+', postfix)
                    if "dipt" in self.hist_name:
                        text = ("{:.0f}".format(float(edges[0])) + r"$\mathit{ < m < }$" +
                                "{:.0f}".format(float(edges[1])) + " GeV")
                    else:
                        text = (r"$\mathit{p_{T} < }$" + "{:.0f}".format(float(edges[1])) + " GeV")
                    data_hist.set_text_for_plot(text)
            return data_hists

    def make_expectation_hist(self, expectation_dict, second_axis_postfix=""):

        total_expectation_hist = None
        for group_name, files in expectation_dict.items():
            # check if hist_prefix needed
            group_name_split = group_name.split(":")
            if len(group_name_split) == 2:
                group_name = group_name_split[0]
                hist_prefix = group_name_split[1]
            elif len(group_name_split) == 1:
                group_name = group_name_split[0]
                hist_prefix = ""
            else:
                # somthing wrong
                return
            group = self.experiment_file_pather.make_mc_group(self.period, self.channel, *files,
                                                              group_name=group_name, hist_prefix=hist_prefix,
                                                              file_path_postfix=self.file_path_postfix,
                                                              hist_path_postfix=self.hist_path_postfix)
            hist = group.get_hist(self.hist_name + second_axis_postfix)
            hist.set_hist_config(self.bin_width_norm, self.axis_steering, self.normalize)

            if total_expectation_hist is None:
                total_expectation_hist = hist
            else:
                total_expectation_hist += hist
        total_expectation_hist.set_hist_config(self.bin_width_norm, self.axis_steering, self.normalize)
        return total_expectation_hist

    def get_total_expectation_hist(self, exp_type='total', force_1d_output=False):

        if exp_type == "total":
            expectation_dict = self.expectation_dict
        elif exp_type == "signal":
            expectation_dict = self.signal_dict
        elif exp_type == 'bg':
            expectation_dict = self.bg_dict
        else:
            return

        if self.second_axis_postfix is None:
            total_expectation_hist = self.make_expectation_hist(expectation_dict)  # dict, hist, second_axis_postfix
            if force_1d_output:
                return total_expectation_hist.get_1d_hists()
            else:
                return total_expectation_hist
        else:
            total_expectation_hists = []
            for postfix in self.second_axis_postfix:
                total_expectation_hist = self.make_expectation_hist(expectation_dict, postfix)
                total_expectation_hists.append(total_expectation_hist)
            return total_expectation_hists
