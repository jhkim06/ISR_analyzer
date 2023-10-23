from Plotter import Plotter
from HistProducer import HistProducer
import os


class Analyzer:

    def __init__(self, file_pather, signal_dict=None, bg_dict=None):

        self.file_pather = file_pather  # enable access to an experiment results
        self.experiment = self.file_pather.experiment_name
        self.signal_dict = signal_dict
        self.bg_dict = bg_dict

        self.period = ''
        self.channel = ''
        self.file_path_postfix = ''
        self.hist_producer = HistProducer(self.file_pather, self.signal_dict, self.bg_dict)
        self.path_for_plots = './plots/' + self.experiment + "/"
        self.plotter = Plotter(self.experiment, self.path_for_plots)

        # make directory
        if not os.path.exists(self.path_for_plots):
            os.mkdir(self.path_for_plots)

    def get_channel_label(self):
        if self.channel == 'ee':
            return 'ee'
        else:
            return r'\mu\mu'

    def set_base_config(self, period, channel,
                        file_path_postfix='',
                        hist_path_postfix=''):

        self.period = period
        self.channel = channel
        self.file_path_postfix = file_path_postfix

        self.hist_producer.set_base_configs(period, channel,
                                            file_path_postfix,
                                            hist_path_postfix)
        self.plotter.set_period_channel(period, channel,
                                        hist_path_postfix)

    def draw_data_expectations(self, hist_name, hist_path_postfix="",
                               bin_width_norm=False, axis_steering="", figsize=(8, 8), **kwargs):

        print(f"{self.period} {self.channel} plotting {hist_name}...")

        # TODO error handle
        if self.period == '' or self.channel == '':
            return
        self.hist_producer.set_configs(hist_name,
                                       hist_path_postfix=hist_path_postfix,
                                       bin_width_norm=bin_width_norm,
                                       axis_steering=axis_steering)

        data_hist = self.hist_producer.get_data_hist()
        total_expectation_hist = self.hist_producer.get_total_expectation_hist()

        self.plotter.make_comparison_plot(total_expectation_hist, data_hist, figsize=figsize)

    # def draw(denominator, *nominators)
    def draw_comparisons(self, denominator, nominator, hist_name, hist_type='data', show_mean_line=False):
        # local HistProducer to get hist
        denominator_producer = HistProducer(self.file_pather, self.signal_dict, self.bg_dict)
        denominator_producer.set_base_configs(denominator[0], denominator[1])
        denominator_producer.set_configs(hist_name,
                                         hist_path_postfix=denominator[2],
                                         bin_width_norm=True,
                                         normalize=True)

        nominator_producer = HistProducer(self.file_pather, self.signal_dict, self.bg_dict)
        nominator_producer.set_base_configs(nominator[0], nominator[1])
        nominator_producer.set_configs(hist_name,
                                       hist_path_postfix=nominator[2],
                                       bin_width_norm=True,
                                       normalize=True)
        if hist_type == 'data':
            denominator_hist = denominator_producer.get_data_hist()
            nominator_hist = nominator_producer.get_data_hist()
        elif hist_type == 'signal':
            denominator_hist = denominator_producer.get_total_expectation_hist(exp_type=hist_type)
            nominator_hist = nominator_producer.get_total_expectation_hist(exp_type=hist_type)

        # draw
        self.plotter.make_comparison_plot(nominator_hist, denominator_hist, show_mean_line=show_mean_line,
                                          ymin=0.8, ymax=1.2)
