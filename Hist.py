import numpy as np
from collections import namedtuple
from collections import OrderedDict
import re
import math


HistNumpy = namedtuple('Hist', ['values', 'bins', 'errors'])
Stack = namedtuple('Stack', ['values_list', 'bins', 'errors_list'])


def create_extracted_hist(hist, axis_steering):
    hist.set_axis_steering(axis_steering)
    raw_hist = hist.hist()  #
    sys_raw_hist_dict = None
    if hist.sys_on:
        sys_raw_hist_dict = dict()
        for sys_name, variations in hist.root_sys_hist_dict.items():
            sys_raw_hist_dict[sys_name] = dict()
            for variation in variations:
                sys_raw_hist_dict[sys_name][variation] = hist.hist(sys_name, variation)
    return raw_hist, sys_raw_hist_dict


def create_projected_hist(hist, axis, first_bin, last_bin, option):
    raw_hist = hist.projected_hist('', '', axis, first_bin, last_bin, option)
    sys_raw_hist_dict = None
    if hist.sys_on:
        sys_raw_hist_dict = dict()
        for sys_name, variations in hist.root_sys_hist_dict.items():
            sys_raw_hist_dict[sys_name] = dict()
            for variation in variations:
                sys_raw_hist_dict[sys_name][variation] = hist.projected_hist(sys_name, variation,
                                                                             axis, first_bin, last_bin, option)
    return raw_hist, sys_raw_hist_dict


class Hist:

    def __init__(self, hist_name, root_hist, root_sys_hist_dict, file_group, multiple_hists=False):

        self.multiple_hists = multiple_hists
        self.hist_name = hist_name
        self.root_hist = root_hist
        self.root_sys_hist_dict = root_sys_hist_dict

        # flag for measurement, expectation
        self.is_mc = True

        self.sys_on = False
        if self.root_sys_hist_dict is not None and len(self.root_sys_hist_dict) > 0:
            self.sys_on = True

        self.unfold_bin = dict()

        self.file_group = file_group
        self.group_label = self.file_group.get_group_name

        self.bin1 = None  #
        self.bin2 = None  #
        if "tunfold" in self.hist_name:
            self.bin1, self.bin2 = self.get_bin_name()
            self.unfold_bin[self.bin1] = self.file_group.get_unfold_bin(self.bin1)
            if self.bin2 is not None:
                self.unfold_bin[self.bin2] = self.file_group.get_unfold_bin(self.bin2)

        self.use_axis_steering = False
        self.axis_steering = ""

        self.bin_width_norm = False
        self.normalize = False

        self.hist_dict = OrderedDict()
        if not self.multiple_hists:
            self.hist_dict[self.group_label] = self

        # text to show in plot
        self.text_for_plot = None

    def set_hist_name(self, hist_name):
        self.hist_name = hist_name

    def set_is_mc(self, is_mc):
        self.is_mc = is_mc

    def update_group_label(self, new_group_label='', postfix=''):
        if postfix == '':
            self.group_label = new_group_label
        else:
            if new_group_label == '':
                self.group_label += postfix
            else:
                self.group_label = new_group_label + postfix

    def set_text_for_plot(self, text):
        self.text_for_plot = text

    def get_text_for_plot(self):
        return self.text_for_plot

    def axis_steering_text(self, text=''):

        if self.unfold_bin[self.bin1].GetDistributionAxisLabel(1) == "dimass":
            edges = self.unfold_bin[self.bin1].GetDistributionBinning(1)
            mass_windows = [edge for edge in edges]

            # get mass index from axis steering
            matches = re.findall(r'(\[[^]]*])', self.axis_steering)
            index = int(re.findall(r'(\d)', matches[1])[0])

            text = ("{:.0f}".format(mass_windows[index]) + r"$\mathit{ < m < }$" +
                    "{:.0f}".format(mass_windows[index + 1]) + " [GeV]")
        elif self.unfold_bin[self.bin1].GetDistributionAxisLabel(1) == "dipt":
            edges = self.unfold_bin[self.bin1].GetDistributionBinning(1)
            pt_cut = edges[1]

            text = (r"$\mathit{p_{T} < }$" +
                    "{:.0f}".format(pt_cut) + " [GeV]")
        else:
            text = text

        return text

    def get_bin_name(self):
        matches = re.findall(r'(\[[^]]*])', self.hist_name)
        var_name = re.sub(r'(\([^]]*\))', '', matches[1])
        bin1 = re.sub(r'\[[^[]*?__', '[', matches[2])
        bin_name1 = "[tunfold-bin]_" + var_name + "_" + bin1
        bin_name2 = None
        if len(matches) > 3:
            bin2 = re.sub(r'\[[^[]*?__', '[', matches[3])
            bin_name2 = "[tunfold-bin]_" + var_name + "_" + bin2
        return bin_name1, bin_name2

    def unset_hist_config(self):
        self.use_axis_steering = False
        self.bin_width_norm = False
        self.normalize = False

    def set_axis_steering(self, axis_steering):

        if axis_steering != "" and self.unfold_bin is not None:
            self.use_axis_steering = True
            self.axis_steering = axis_steering

            self.text_for_plot = self.axis_steering_text()
        else:
            self.use_axis_steering = False

    def set_hist_config(self, bin_width_norm, axis_steering, normalize=False):  # fix histogram configuration
        self.set_axis_steering(axis_steering)
        self.bin_width_norm = bin_width_norm
        self.normalize = normalize

        if self.multiple_hists:
            for _, hist in self.hist_dict.items():
                hist.set_hist_config(bin_width_norm, axis_steering, normalize)

    def get_unfold_bin(self, bin_name):
        return self.unfold_bin[bin_name]

    def get_root_hist(self, sys_name="", sys_variation=""):
        if sys_name == "" and sys_variation == "":
            return self.root_hist
        else:
            if self.root_sys_hist_dict is not None:
                if sys_name in self.root_sys_hist_dict:
                    if sys_variation in self.root_sys_hist_dict[sys_name]:
                        return self.root_sys_hist_dict[sys_name][sys_variation]
                    else:
                        return self.root_hist
                else:
                    return self.root_hist
            else:
                return self.root_hist

    def get_sqrt_sys_hist(self):
        pass

    # get 1D hists from a 2D tunfold hist
    def get_1d_hists(self):

        out_hists = []
        base_steering = {'dipt': 'O', 'dimass': 'UO'}  # note fixed axis steering
        axis1_name = str(self.unfold_bin[self.bin1].GetDistributionAxisLabel(0))
        axis2_name = str(self.unfold_bin[self.bin1].GetDistributionAxisLabel(1))

        edges = self.unfold_bin[self.bin1].GetDistributionBinning(1)
        second_axis_edges = [edge for edge in edges]

        self_axis_steering = self.axis_steering
        for index in range(len(second_axis_edges) - 1):
            axis_steering = (f"{axis1_name}[{base_steering[axis1_name]}];" +
                             f"{axis2_name}[{base_steering[axis2_name]}C{index}]")

            # create extracted Hist/
            raw_hist, sys_raw_hist_dict = create_extracted_hist(self, axis_steering)
            out_hist = Hist("extracted_" + self.hist_name + "_" + axis1_name + "_" + str(index),
                            raw_hist, sys_raw_hist_dict, self.file_group)
            out_hist.is_mc = self.is_mc
            out_hist.update_group_label(self.group_label)  #
            out_hist.text_for_plot = self.text_for_plot
            if self.multiple_hists:
                extracted_hist_dict = dict()
                for label, hist in self.hist_dict.items():
                    temp_hist, temp_sys_hist_dict = create_extracted_hist(hist, axis_steering)
                    extracted_hist_dict[label] = Hist("extracted_" + axis1_name, temp_hist, temp_sys_hist_dict,
                                                      hist.file_group)
                    extracted_hist_dict[label].is_mc = hist.is_mc
                    extracted_hist_dict[label].update_group_label(self.group_label)
                    hist.set_axis_steering(self_axis_steering)
                out_hist.multiple_hists = True
                out_hist.hist_dict = extracted_hist_dict

            out_hists.append(out_hist)
        self.set_axis_steering(self_axis_steering)
        return out_hists

    # for TH2D
    def get_projected_hist(self, axis='x', first_bin=0, last_bin=-1, option='e'):
        raw_hist, sys_raw_hist_dict = create_projected_hist(self, axis, first_bin, last_bin, option)
        matches = re.findall(r'(\[[^]]*])', self.hist_name)
        new_hist_name = self.hist_name
        if len(matches) == 4:
            new_hist_name = matches[0] + "_" + matches[1] + "_"
            if axis == 'x':
                new_hist_name += matches[3]
            else:
                new_hist_name += matches[2]
        out_hist = Hist("projected_" + axis + "_" + new_hist_name,
                        raw_hist, sys_raw_hist_dict, self.file_group)

        out_hist.is_mc = self.is_mc
        out_hist.update_group_label(self.group_label)  #
        out_hist.text_for_plot = self.text_for_plot
        # TODO multiple hists?

        return out_hist

    # squared root sum for systematics
    def __add__(self, other):

        root_hist = self.root_hist.Clone("clone")
        root_hist.Add(other.root_hist)

        root_sys_hist_dict = self.__make_sys_hist_dict__(other)

        new_hist = Hist(self.hist_name, root_hist, root_sys_hist_dict, self.file_group, True)
        new_hist.set_hist_config(False, self.axis_steering, False)  # set bin_width_norm False

        new_hist.hist_dict.update(other.hist_dict)
        new_hist.hist_dict.update(self.hist_dict)

        return new_hist

    def __loop_sys_dict__(self, output_dict, sys_names, sys_hist_dict, other):

        for sys_name in sys_names:
            output_dict[sys_name] = dict()
            for variation in sys_hist_dict[sys_name]:

                output_dict[sys_name][variation] = self.get_root_hist(sys_name, variation).Clone('clone')
                other_sys_hist = other.get_root_hist(sys_name, variation)

                output_dict[sys_name][variation].Add(other_sys_hist)

    def __make_sys_hist_dict__(self, other):

        if self.sys_on and other.sys_on:
            root_sys_hist_dict = dict()

            common_systematic_names = set(self.root_sys_hist_dict.keys()).intersection(
                set(other.root_sys_hist_dict.keys()))
            self_only_systematic = set(self.root_sys_hist_dict.keys()).difference(
                set(other.root_sys_hist_dict.keys()))
            other_only_systematic = set(other.root_sys_hist_dict.keys()).difference(
                set(self.root_sys_hist_dict.keys()))

            # common systematics
            self.__loop_sys_dict__(root_sys_hist_dict, common_systematic_names, self.root_sys_hist_dict, other)
            # self
            self.__loop_sys_dict__(root_sys_hist_dict, self_only_systematic, self.root_sys_hist_dict, other)
            # other
            self.__loop_sys_dict__(root_sys_hist_dict, other_only_systematic, other.root_sys_hist_dict, other)

        elif self.sys_on and not other.sys_on:
            root_sys_hist_dict = dict()
            self.__loop_sys_dict__(root_sys_hist_dict, self.root_sys_hist_dict.keys(), self.root_sys_hist_dict, other)

        elif not self.sys_on and other.sys_on:
            root_sys_hist_dict = dict()
            self.__loop_sys_dict__(root_sys_hist_dict, other.root_sys_hist_dict.keys(), other.root_sys_hist_dict, other)
        else:
            return None

        return root_sys_hist_dict

    def multiply(self, other):

        # use raw histogram
        root_hist = self.root_hist.Clone("ratio")
        if len(root_hist.GetSumw2()) == 0:
            root_hist.Sumw2()
        root_hist.Multiply(other.root_hist)

        if self.sys_on:
            root_sys_hist_dict = dict()
            for sys_name, variations in self.root_sys_hist_dict.items():
                root_sys_hist_dict[sys_name] = dict()
                for variation in variations:
                    root_sys_hist_dict[sys_name][variation] = \
                        self.root_sys_hist_dict[sys_name][variation].Clone("clone")
                    root_sys_hist_dict[sys_name][variation].Multiply(other.root_hist)
        else:
            root_sys_hist_dict = None

        new_hist = Hist(self.hist_name, root_hist, root_sys_hist_dict, self.file_group)
        new_hist.set_hist_config(False, self.axis_steering, False)  # set bin_width_norm False
        new_hist.is_mc = self.is_mc
        return new_hist

    def divide(self, other=None):

        if other is None:
            denominator = self
        else:
            denominator = other

        # use raw histogram
        root_hist = self.hist()
        denominator_hist = denominator.hist()
        if len(root_hist.GetSumw2()) == 0:
            root_hist.Sumw2()
        root_hist.Divide(denominator_hist)

        if self.sys_on:
            root_sys_hist_dict = dict()
            for sys_name, variations in self.root_sys_hist_dict.items():
                root_sys_hist_dict[sys_name] = dict()
                for variation in variations:
                    # root_sys_hist_dict[sys_name][variation] = \
                    #    self.root_sys_hist_dict[sys_name][variation].Clone("clone")
                    root_sys_hist_dict[sys_name][variation] = self.hist(sys_name, variation)
                    root_sys_hist_dict[sys_name][variation].Divide(denominator_hist)
        else:
            root_sys_hist_dict = None

        new_hist = Hist(self.hist_name, root_hist, root_sys_hist_dict, self.file_group)
        new_hist.set_hist_config(False, self.axis_steering, False)  # set bin_width_norm False
        new_hist.is_mc = self.is_mc
        return new_hist

    def hist(self, sys_name="", sys_variation=""):

        if sys_name == "":
            hist = self.root_hist.Clone("clone")
        else:
            hist = self.root_sys_hist_dict[sys_name][sys_variation].Clone("clone")

        if self.use_axis_steering:
            hist = self.unfold_bin[self.bin1].ExtractHistogram("extracted",
                                                               hist, 0, True, self.axis_steering)
        if self.bin_width_norm:
            if self.normalize:
                hist.Scale(1/hist.Integral(), "width")
            else:
                hist.Scale(1, "width")
        return hist

    def projected_hist(self, sys_name="", sys_variation="",
                       axis='x', first_bin=0, last_bin=-1, option='e'):

        self_use_axis_steering = self.use_axis_steering
        self_bin_width_norm = self.bin_width_norm
        self.use_axis_steering = False
        self.bin_width_norm = False
        matrix = self.hist(sys_name, sys_variation)

        if axis == 'x':
            projected_hist = matrix.ProjectionX('projected_x', first_bin, last_bin, option)
        elif axis == 'y':
            projected_hist = matrix.ProjectionY('projected_y', first_bin, last_bin, option)
        else:
            print("use valid axis, currently {axis} not available.")
            return

        self.use_axis_steering = self_use_axis_steering
        self.bin_width_norm = self_bin_width_norm
        if self.use_axis_steering:
            projected_hist = self.unfold_bin[self.bin1].ExtractHistogram("extracted",
                                                                         projected_hist, 0, True, self.axis_steering)
        if self.bin_width_norm:
            projected_hist.Scale(1, "width")

        return projected_hist

    def get_mean(self, x_min=0, x_max=0, binned_mean=False):
        bin_width_norm = self.bin_width_norm
        self.bin_width_norm = False

        hist = self.hist()
        if binned_mean:
            min_axis = hist.GetXaxis().GetXmin()
            max_axis = hist.GetXaxis().GetXmax()
            hist.GetXaxis().SetRangeUser(min_axis, max_axis)
        if x_min < x_max:
            hist.GetXaxis().SetRangeUser(x_min, x_max)

        mean = hist.GetMean()

        errors = dict()
        errors["stat"] = hist.GetMeanError()

        # loop over sys names, and get symmetric error
        # FIXME update to use appropriate method for each systematic
        if self.sys_on:
            for sys_name, variations in self.root_sys_hist_dict.items():
                error_list = []
                for variation in variations:
                    sys_hist = self.hist(sys_name, variation)
                    if x_min < x_max:
                        sys_hist.GetXaxis().SetRangeUser(x_min, x_max)
                    sys_mean = sys_hist.GetMean()

                    error_list.append(sys_mean-mean)
                errors[sys_name] = math.sqrt(np.sum(np.square(error_list)))

        self.bin_width_norm = bin_width_norm
        hist.GetXaxis().SetRangeUser(0, 0)  # unset the range
        return mean, errors

    # return symmetry systematic in numpy format
    def get_sys_numpy(self, sys_name="total"):
        default, _, _ = self.to_numpy()

        if sys_name == "total":
            # take absolute average variation for each systematic and square root sum
            # if systematic only have one variation don't use it
            total_systematic = None
            for sys_name, variations in self.root_sys_hist_dict.items():
                systematic_variations = None

                for variation in variations:
                    systematic, _, _ = self.to_numpy(sys_name, variation)
                    systematic = systematic - default
                    systematic = np.expand_dims(np.absolute(systematic), axis=0)
                    # update systematic variation
                    if systematic_variations is None:
                        systematic_variations = systematic
                    else:
                        systematic_variations = np.append(systematic_variations, systematic, axis=0)

                # get mean
                systematic_variations = np.mean(systematic_variations, axis=0)
                systematic_variations = np.expand_dims(systematic_variations, axis=0)

                if total_systematic is None:
                    total_systematic = systematic_variations
                else:
                    total_systematic = np.append(total_systematic, systematic_variations, axis=0)

            total_systematic = np.sqrt(np.sum(np.square(total_systematic), axis=0))
            return total_systematic

        elif sys_name == "systematic":
            pass
        elif sys_name in self.root_sys_hist_dict:
            pass
        else:
            pass

    # THistogram to numpy
    def to_numpy(self, sys_name="", sys_variation=""):

        values = []
        bins = []
        error = []

        hist = self.hist(sys_name, sys_variation)
        # check if bins have labels
        bin_labels = hist.GetXaxis().GetLabels()
        n_bins_x = hist.GetNbinsX()
        for i_bin in range(n_bins_x):
            values.append(hist.GetBinContent(i_bin + 1))
            error.append(hist.GetBinError(i_bin + 1))
            if bin_labels:
                bins.append(str(bin_labels[i_bin]))
            else:
                bins.append(hist.GetXaxis().GetBinLowEdge(i_bin + 1))

        if bin_labels:
            pass
        else:
            bins.append(bins[-1] + hist.GetBinWidth(n_bins_x))

        values = np.array(values)
        bins = np.array(bins)
        errors = np.array(error)  # stat

        h = HistNumpy(values, bins, errors)
        return h

    def get_stack(self):

        values_list = []
        bins = None
        errors_list = []
        labels = []

        for label, hist in self.hist_dict.items():
            hist_np = hist.to_numpy()
            values_list.append(hist_np.values)
            if bins is None:
                bins = hist_np.bins
            errors_list.append(hist_np.errors)
            labels.append(label)

        stack = Stack(values_list, bins, errors_list)
        return stack, labels

    # return dictionary of numpy
    def to_sys_numpy(self):  # to_sys_numpy_dict()

        sys_numpy_dict = dict()
        for sys_name, variations in self.root_sys_hist_dict.items():
            sys_numpy_dict[sys_name] = dict()
            for variation in variations:
                sys_numpy_dict[sys_name][variation] = self.to_numpy(sys_name, variation)

        return sys_numpy_dict
