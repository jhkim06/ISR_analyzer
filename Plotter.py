import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import (FixedLocator, FixedFormatter)
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import numpy as np
from matplotlib.offsetbox import AnchoredText
import os


def adjust_x_lim(axis, bins):
    if bins[-1] >= 1e3:
        axis.set_xscale("log")
        if bins[0] == 0:
            axis.set_xlim(0.2, bins[-1])
        else:
            axis.set_xlim(bins[0], bins[-1])

            minor_x_ticks = FixedLocator([50, 60, 70, 80, 90, 200, 300, 400, 500, 600, 700, 800, 900])
            minor_x_ticklabels = ['50', '60', '70', '', '90', '200', '', '', '', '', '', '', '']
            axis.xaxis.set_minor_locator(minor_x_ticks)
            axis.xaxis.set_minor_formatter(FixedFormatter(minor_x_ticklabels))

            major_x_ticks = FixedLocator([100, 1000])
            major_x_ticklabels = ['', '1000']
            axis.xaxis.set_major_locator(major_x_ticks)
            axis.xaxis.set_major_formatter(FixedFormatter(major_x_ticklabels))
            axis.tick_params(axis='x', which='both', labelsize=30)
    else:
        axis.set_xlim(bins[0], bins[-1])


def make_error_boxes(default_value, bins, errors):  # TODO add options for hatch cosmetics

    # artist[0][0][0].get_data()
    # artist[0][0][2][0].get_segments()  # bin content
    # artist[0][0][2][1].get_segments()  # error boundary

    center = bins[:-1] + np.diff(bins) / 2.
    x_width = np.expand_dims(np.diff(bins) / 2., axis=0)
    x_width = np.append(x_width, x_width, axis=0)

    errors = np.expand_dims(errors, axis=0)
    errors = np.append(errors, errors, axis=0)

    error_boxes = [Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(),
                             fill=False, edgecolor=None, facecolor=None)
                   for x, y, xe, ye in zip(center, default_value, x_width.T, errors.T)]
    pc = PatchCollection(error_boxes, match_original=True, hatch="///", linewidth=0.0, zorder=100)
    # for legend only
    legend_ = Rectangle((0, 0), 0.0, 0.0,
                        fill=False, edgecolor=None, facecolor=None, label='Sys.', hatch="///", linewidth=0.5)
    return pc, legend_


class Plotter:
    def __init__(self, experiment, base_output_dir, **kwargs):

        # mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["royalblue", "cornflowerblue", "turquoise",
        #                                                     'tomato', "salmon"])
        self.experiment = experiment
        self.period = ''
        self.channel = ''
        hep.style.use(self.experiment)
        plt.ioff()  # interactive mode off; not to show figures on creation
        plt.rcParams.update({
            "text.usetex": True,
            # "font.family": "Computer Modern Roman",
            'text.latex.preamble': r'\usepackage{newtxmath}'
        })
        plt.rcParams['axes.linewidth'] = 2.0
        plt.rcParams['hatch.linewidth'] = 0.5

        self.rows = 1
        self.cols = 1

        self.fig = None
        self.axs = None
        self.current_axis = None

        self.hists = []
        self.hist_kwargs = []
        self.hist_draw_option = []
        self.denominator_index = -1

        self.text = ''
        self.use_text_as_label = False
        self.text_written = False
        self.use_x_axis_tick_labels_from_hist = False

        # make directory to save plots
        self.base_output_dir = base_output_dir
        self.out_dir = ''

    def reset(self):

        self.rows = 1
        self.cols = 1

        self.fig = None
        self.axs = None
        self.current_axis = None

        self.hists.clear()
        self.hist_kwargs.clear()
        self.hist_draw_option.clear()
        self.denominator_index = -1

        self.text = ''
        self.use_text_as_label = False
        self.text_written = False
        self.use_x_axis_tick_labels_from_hist = False

    def set_period_channel(self, period, channel, hist_path_postfix=''):
        self.period = period
        self.channel = channel

        self.out_dir = self.base_output_dir + self.period + "/" + self.channel + "/"
        if hist_path_postfix != '':
            self.out_dir += hist_path_postfix + "/"
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

    def set_text(self, raw_hist=None, text=None):
        # TODO allow another text
        if raw_hist:
            self.text = raw_hist.get_text_for_plot()
        if text:
            self.text = text

    def create_subplots(self, rows=1, cols=1, figsize=(8, 8), **gridspec_kw):
        if rows != self.rows:
            self.rows = rows
        if cols != self.cols:
            self.cols = cols

        self.fig, self.axs = plt.subplots(self.rows, self.cols,
                                          figsize=figsize,
                                          gridspec_kw=gridspec_kw)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, hspace=0.05)

    def add_hist(self, hist, is_denominator=False, as_stack=False, **kwargs):
        self.hists.append(hist)
        self.hist_kwargs.append(kwargs)
        self.hist_draw_option.append(as_stack)
        if is_denominator:
            self.denominator_index = len(self.hists) - 1

    def set_current_axis(self, row, col):
        if self.rows == 1 and self.cols == 1:
            self.current_axis = self.axs
            return self.axs
        elif self.rows == 1 or self.cols == 1:
            if self.rows == 1:
                index = col
            else:
                index = row
            self.current_axis = self.axs[index]
            return self.axs[index]
        else:
            self.current_axis = self.axs[row][col]
            return self.axs[row][col]

    def get_axis(self, row, col):
        return self.set_current_axis(row, col)

    def set_y_axis_label(self, hist=None, label=None):
        axis = self.current_axis
        fontsize = 30
        if label is None:  # determine based on Hist
            if hist.normalize:
                axis.set_ylabel("a.u.", fontsize=fontsize)
            elif hist.bin_width_norm:
                axis.set_ylabel("Events/1 GeV", fontsize=fontsize)
            else:
                axis.set_ylabel("Events/bin", fontsize=fontsize)
        else:
            axis.set_ylabel(label, fontsize=fontsize)

    def set_x_axis_label(self, hist=None, label=None):
        axis = self.current_axis
        fontsize = 30
        if label is None:
            axis.set_xlabel(hist.hist_name, fontsize=fontsize)
        else:
            axis.set_xlabel(label, fontsize=fontsize)

    def x_axis_tick_labels(self, hist, bins):
        if self.use_x_axis_tick_labels_from_hist:
            adjust_x_lim(self.current_axis, bins)
            self.set_tick_labels(bins, hist.bins)
        else:
            adjust_x_lim(self.current_axis, hist.bins)

    def save_fig(self, out_name=''):
        if out_name == '':
            out_file_name = self.hists[0].hist_name
        else:
            out_file_name = out_name

        if self.text:
            out_file_name += "_" + self.text
        print(f"save plot... {out_file_name}")
        self.save_plot(self.out_dir + out_file_name)
        self.reset()
        plt.close()

    def draw(self):
        for index, hist in enumerate(self.hists):
            if hist.multiple_hists and self.hist_draw_option[index]:
                self.draw_stack(hist, use_mplhep=False)
            else:
                artist = self.draw_hist(hist, **self.hist_kwargs[index])
                if "histtype" in self.hist_kwargs[index]:
                    if self.hist_kwargs[index]["histtype"] == "step":
                        ec = artist[0][0].get_edgecolor()
                        self.hist_kwargs[index]['color'] = ec

                    elif self.hist_kwargs[index]["histtype"] == "errorbar":
                        mec = artist[0][0][0].get_markeredgecolor()
                        self.hist_kwargs[index]["mec"] = mec
                        mfc = artist[0][0][0].get_markerfacecolor()
                        self.hist_kwargs[index]["mfc"] = mfc

                    else:  # self.hist_kwargs[index]["histtype"] == "fill":
                        ec = artist[0][0].get_edgecolor()
                        self.hist_kwargs[index]["ec"] = ec
                        fc = artist[0][0].get_facecolor()
                        self.hist_kwargs[index]["fc"] = fc
                else:
                    ec = artist[0][0].get_edgecolor()
                    self.hist_kwargs[index]["ec"] = ec

            if index == 0:
                self.set_y_axis_label(hist)
                self.set_x_axis_label(hist)

    def draw_ratio(self, y_axis_label=None):
        denominator_hist = self.hists[self.denominator_index]
        is_first = True
        denominator_colr = 'black'
        for index, hist in enumerate(self.hists):
            if index == self.denominator_index:
                if 'color' in self.hist_kwargs[index]:
                    denominator_colr = self.hist_kwargs[index]['color']
                elif "mec" in self.hist_kwargs[index]:
                    denominator_colr = self.hist_kwargs[index]['mec']
                elif "mfc" in self.hist_kwargs[index]:
                    denominator_colr = self.hist_kwargs[index]['mfc']
                elif "ec" in self.hist_kwargs[index]:
                    denominator_colr = self.hist_kwargs[index]['ec']
                elif "fc" in self.hist_kwargs[index]:
                    denominator_colr = self.hist_kwargs[index]['fc']
                else:
                    pass
                continue
            ratio = hist.divide(denominator_hist)
            self.draw_hist(ratio, **self.hist_kwargs[index])
            if is_first:
                label = ratio.hist_name
                if y_axis_label:
                    label = y_axis_label
                self.set_y_axis_label(ratio, label)
                is_first = False
        self.current_axis.axhline(y=1, color=denominator_colr, linestyle='--', linewidth=1)

    def set_experiment_label(self, label="Preliminary"):
        self.current_axis.draw(self.current_axis.figure.canvas.get_renderer())

        plt.rcParams['text.usetex'] = False
        hep.cms.label(label, data=True, year=self.period, ax=self.current_axis, loc=0, pad=.0)
        plt.rcParams['text.usetex'] = True

    def set_axis_config(self, do_mpl_magic=True, set_legend=True, legend_loc='best', sort_legend=True):
        if set_legend:
            # check n hist
            if len(self.hists) > 10:
                ncol = 2
                fontsize = 20
            else:
                ncol = 1
                fontsize = 25
            hep.plot.hist_legend(self.current_axis, loc=legend_loc,
                                 ncol=ncol, fontsize=fontsize)
            if sort_legend:
                handles, labels = self.current_axis.get_legend_handles_labels()
                hep.plot.sort_legend(ax=self.current_axis, order=labels[::-1])

        if do_mpl_magic:
            hep.plot.yscale_legend(self.current_axis)
            if self.text_written:
                # hep.plot.yscale_text(self.current_axis)
                hep.plot.mpl_magic(self.current_axis)

    def set_y_axis_config(self, label=None, log_scale=False,
                          remove_first_tick_label=False,
                          set_min=None, set_max=None,
                          set_grid=False):
        self.current_axis.tick_params(axis='y', labelsize=30)
        if label:
            self.set_y_axis_label(label=label)

        if log_scale:
            self.current_axis.set_yscale("log")
        else:
            self.current_axis.ticklabel_format(style='sci', axis='y', scilimits=(0, 4))

        if set_min is not None:
            self.current_axis.set_ylim(ymin=set_min)

        if set_max is not None:
            self.current_axis.set_ylim(ymax=set_max)

        if remove_first_tick_label:
            self.current_axis.yaxis.get_major_ticks()[0].set_visible(False)

        if set_grid:
            self.current_axis.grid(axis='y', which='both')

    def set_x_axis_config(self, remove_tick_labels=False, set_grid=False, label=None, log_scale=False):

        self.current_axis.tick_params(axis='x', which='both', labelsize=30)
        if label:
            self.set_x_axis_label(label=label)

        if remove_tick_labels:
            self.current_axis.set_xticklabels([])  # delete major labels
            self.current_axis.set_xticklabels([], minor=True)  # delete minor labels
            self.set_x_axis_label(label="")

        if set_grid:
            self.current_axis.grid(axis='x', which='both')

        if log_scale:
            x_min, x_max = self.current_axis.get_xlim()
            if x_min == 0:
                self.current_axis.set_xlim(0.2, x_max)
            self.current_axis.set_xscale("log")

    # FIXME move this method to Hist
    def get_x_bins(self, hist):
        if all(isinstance(n, str) for n in hist.bins):
            self.use_x_axis_tick_labels_from_hist = True
            bins = list(range(len(hist.bins) + 1))
        else:
            bins = hist.bins
        return bins

    def draw_hist(self, raw_hist, **kwargs):
        hist = raw_hist.to_numpy()
        bins = self.get_x_bins(hist)

        label = raw_hist.hist_label
        if self.use_text_as_label:
            label = raw_hist.get_text_for_plot()
        else:
            if self.text == "":  # update only when self.text is not set
                self.set_text(raw_hist=raw_hist)

        if 'label' in kwargs:
            label = kwargs['label']
            kwargs.pop('label')

        artist = hep.histplot((hist.values, bins), ax=self.current_axis, yerr=hist.errors,
                              label=label, **kwargs)
        self.x_axis_tick_labels(hist, bins)

        if self.text != '' and not self.text_written:
            at = AnchoredText(self.text, loc='upper left', prop=dict(size=20), frameon=False)
            at.patch.set_boxstyle("round, pad=0., rounding_size=0.2")
            self.current_axis.add_artist(at)
            self.text_written = True

        if raw_hist.sys_on and raw_hist.is_mc:  # draw error band for total uncertainty
            sys_error = raw_hist.get_sys_numpy()
            pc, legend_ = make_error_boxes(hist.values, hist.bins, sys_error)
            # TODO option to choose error boxes and simple error bar
            self.current_axis.add_collection(pc)
            self.current_axis.add_patch(legend_)

        return artist

    def draw_stack(self, raw_hist, use_mplhep=True, **kwargs):
        stack, labels = raw_hist.get_stack()
        bins = self.get_x_bins(stack)

        if use_mplhep:
            hep.histplot(stack.values_list, ax=self.current_axis, stack=True, histtype="fill", bins=bins,
                         **kwargs)
        else:
            bottom = 0
            for index, values in enumerate(stack.values_list):
                self.current_axis.hist(bins[:-1], bins, weights=values, histtype='bar',
                                       bottom=bottom, label=labels[index])
                bottom += values
        self.x_axis_tick_labels(stack, bins)

        if raw_hist.sys_on:
            hist = raw_hist.to_numpy()
            sys_error = raw_hist.get_sys_numpy()
            pc, legend_ = make_error_boxes(hist.values, hist.bins, sys_error)
            self.current_axis.add_collection(pc)
            self.current_axis.add_patch(legend_)

    def set_tick_labels(self, bins, labels):
        self.current_axis.minorticks_off()
        major_ticks = FixedLocator(bins[1:])
        self.current_axis.xaxis.set_major_locator(major_ticks)
        self.current_axis.xaxis.set_major_formatter(FixedFormatter(labels))
        self.current_axis.tick_params(axis='x', labelrotation=45)

    # FIXME draw_graph()
    def draw_isr_data_frame(self, dimass_df, dipt_df,
                            ymin=14, ymax=30, xmin=55, xmax=300,
                            new_fig=True, **kwargs):
        if new_fig:
            self.rows = 1
            self.cols = 1

            self.fig, self.axs = plt.subplots(self.rows, self.cols, figsize=(11, 8))
            self.set_current_axis(0, 0)
            plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9, hspace=0.05)
            self.set_isr_plot_cosmetics(ymin, ymax, xmin, xmax)

        stat_type = "mean"
        error_name = "total error"

        self.set_experiment_label(label="Preliminary")
        self.current_axis.errorbar(dimass_df[stat_type], dipt_df[stat_type],
                                   xerr=dimass_df[error_name], yerr=dipt_df[error_name],
                                   **kwargs)
        hep.plot.hist_legend(self.current_axis, loc='best')

        final_plot_name = self.out_dir + "isr_result"
        print(f"save plot... {final_plot_name}")
        self.save_plot(final_plot_name)  # fixme save path
        plt.close()

    def set_isr_plot_cosmetics(self, ymin=14, ymax=30, xmin=55, xmax=300):

        axis = self.current_axis
        axis.set_xscale("log")
        axis.set_ylim(ymin, ymax)
        axis.set_xlim(xmin, xmax)

        minor_x_ticks = FixedLocator([50, 60, 70, 80, 90, 200])
        minor_x_ticklabels = ['50', '60', '70', '80', '90', '200']
        axis.xaxis.set_minor_locator(minor_x_ticks)
        axis.xaxis.set_minor_formatter(FixedFormatter(minor_x_ticklabels))

        major_x_ticks = FixedLocator([100])
        major_x_ticklabels = ['']
        axis.xaxis.set_major_locator(major_x_ticks)
        axis.xaxis.set_major_formatter(FixedFormatter(major_x_ticklabels))

        axis.set_ylabel("Mean transverse momentum [GeV]", fontsize=20)
        axis.set_xlabel("Mean mass [GeV]", fontsize=20)

        axis.grid(axis='x', which='both')
        axis.grid(axis='y')

    def draw_isr_error_df(self, df):
        self.rows = 1
        self.cols = 1

        self.fig, self.axs = plt.subplots(self.rows, self.cols, figsize=(15, 8))
        plt.subplots_adjust(left=0.15, right=0.7, bottom=0.15, top=0.9, hspace=0.05)
        axis = self.set_current_axis(0, 0)
        hep.cms.label("Preliminary", data=True, year=self.period, ax=axis, fontsize=20, loc=0)

        df.iloc[:, 2:].plot(kind='bar', ax=self.current_axis)
        handles, legend_labels = axis.get_legend_handles_labels()
        axis.legend(handles, legend_labels, bbox_to_anchor=(1.1, 1.0), loc='upper left', borderaxespad=0.,
                    fontsize=18)

        # show only top 3
        top_n_error = 3
        n_label = len(axis.containers)
        n_mass_index = len(axis.containers[0])
        max_values = []
        for mass_index in range(n_mass_index):
            # print(f'{mass_index} mass window')
            temp_list = []
            for label_index in range(n_label):
                temp_list.append(axis.containers[label_index][mass_index].get_height())
            temp_np = np.array(temp_list)
            index_list = np.argsort(temp_np)
            if len(index_list) > top_n_error:
                max_values.append((temp_np[index_list[-(top_n_error+1)]], temp_np[index_list[-1]]))
            else:
                max_values.append((temp_np[index_list[-1]], temp_np[index_list[-1]]))

        # loop for mass index
        for c, l in zip(axis.containers, legend_labels):
            # customize the labels: only plot values greater than 0 and append the legend label
            labels = [
                f'{l}' if (v.get_height()) > max_values[index][0] and (v.get_height()) != max_values[index][1] else ''
                for index, v in enumerate(c)]  # loop over mass index
            axis.bar_label(c, labels=labels, label_type='edge', rotation=90, fontsize=10)

        axis.set_y_axis_label("mean transverse momentum Unc. [GeV]", fontsize=20)
        axis.set_x_axis_label("mass window index")
        axis.tick_params(axis='x', labelrotation=0)
        axis.minorticks_off()

        final_plot_name = self.out_dir + "isr_error"
        print(f"save plot... {final_plot_name}")
        self.save_plot(final_plot_name)  # fixme save path
        plt.close()

    def clear_all_axis(self):

        plt.rcdefaults()
        hep.style.use(self.experiment)
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['hatch.linewidth'] = 0.5

        if self.rows == 1 and self.cols == 1:
            self.axs.cla()
        elif self.rows == 1 or self.cols == 1:
            if self.rows == 1:
                length = self.cols
            else:
                length = self.rows
            for index in range(length):
                self.axs[index].cla()
        else:
            for row in range(self.rows):
                for col in range(self.cols):
                    self.axs[row][col].cla()

    def delete_axs(self):
        pass

    def save_plot(self, plot_name):
        self.fig.savefig(plot_name+".pdf")
        plt.close()
