from TObjectPather import TObjectPather
import ROOT
import uproot
from Hist import Hist


class FileGroup:
    def __init__(self, exp_name, period_name, channel_name, group_name, file_type,  # FIXME check if file_type used
                 file_dict, sys_file_dict=None,
                 hist_path_postfix="", hist_prefix=""):

        self.experiment_name = exp_name
        self.file_type = file_type
        self.period_name = period_name
        self.channel_name = channel_name
        self.hist_path_postfix = hist_path_postfix

        self.group_name = group_name
        self.is_mc = True

        self.hist_prefix = ""
        self.file_dict = file_dict  # dictionary of file paths
        self.sys_file_dict = sys_file_dict

        # pather for objects in files (nominal histogram)
        self.object_pather = TObjectPather(period_name, channel_name,
                                           hist_path_postfix, hist_prefix)
        self.sys_available = False
        if sys_file_dict is not None:
            self.sys_available = True

    @property
    def get_group_name(self):
        return self.group_name

    def set_is_mc(self, is_mc):
        self.is_mc = is_mc

    def get_sys_hist_postfix(self, hist_name):
        if not self.sys_available:
            return None
        else:
            sys_file_path = self.sys_file_dict[list(self.sys_file_dict.keys())[0]]
            file = uproot.open(sys_file_path)
            # assumed there is ';' end of line
            sys_hist_postfix = []
            for key in file.keys():
                if hist_name in key:
                    sys_postfix = key.split(hist_name)[1].split(";")[0]
                    if sys_postfix != '':
                        sys_hist_postfix.append(sys_postfix)

            if len(sys_hist_postfix) == 0:
                return None
            else:
                return sys_hist_postfix

    def add_hists(self, hist_path, is_sys=False):
        hist_total = None
        ROOT.TH1.AddDirectory(False)
        # loop over files and open histogram with hist path and add them
        if is_sys:
            file_dict = self.sys_file_dict
        else:
            file_dict = self.file_dict

        for file_label, file_path in file_dict.items():
            file = ROOT.TFile.Open(file_path, 'r')
            hist = file.Get(hist_path)
            file.Close()

            if hist_total is None:
                hist_total = hist.Clone(self.group_name)
            else:
                hist_total.Add(hist)

        return hist_total

    def get_unfold_bin(self, bin_name):

        # bin_path = self.object_pather.get_path() + bin_name  # FIXME
        bin_path = bin_name

        file_path = self.file_dict[list(self.file_dict.keys())[0]]
        file = ROOT.TFile.Open(file_path, 'r')
        unfold_bin = file.Get(bin_path)
        file.Close()

        return unfold_bin

    def get_hist(self, hist_name):
        hist_path = self.object_pather.get_hist_path() + hist_name

        ROOT.TH1.AddDirectory(False)
        hist_total = self.add_hists(hist_path)

        # get list of systematic labels
        sys_hist_dict = None
        sys_hist_postfixs = self.get_sys_hist_postfix(hist_path)

        if sys_hist_postfixs is not None:
            sys_names = list({prefix.split("_")[1] for prefix in sys_hist_postfixs})
            sys_hist_dict = dict()
            for sys_name in sys_names:
                sys_hist_dict[sys_name] = dict()

            for hist_postfix in sys_hist_postfixs:
                sys_name = hist_postfix.split("_")[1]
                sys_variation = "_".join(hist_postfix.split("_")[2:])

                sys_hist_path = self.object_pather.get_hist_path(sys=True) + hist_name + hist_postfix
                sys_hist_total = self.add_hists(sys_hist_path, is_sys=True)

                sys_hist_dict[sys_name][sys_variation] = sys_hist_total

        hist = Hist(hist_name, hist_total, sys_hist_dict, self)
        if self.is_mc is False:
            hist.set_is_mc(False)

        return hist
