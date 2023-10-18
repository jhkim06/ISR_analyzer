import json
from FileGroup import FileGroup


class FilePather:
    def __init__(self, base_dir, data_config_json, mc_config_json, exp_name="CMS", sys_on=False):
        self.experiment_name = exp_name
        self.base_dir = base_dir
        self.input_files = dict()

        # file paths
        with open(self.base_dir + data_config_json, 'r') as config_file:
            config_json = json.load(config_file)
            self.input_files["data"] = config_json
        with open(self.base_dir + mc_config_json, 'r') as config_file:
            config_json = json.load(config_file)
            self.input_files["mc"] = config_json
        self.sys_on = sys_on

    def get_data_path(self, period_name, channel_name, file_label, file_path_postfix=""):
        file_name = self.input_files["data"][channel_name][period_name][file_label]
        path = self.base_dir + "/".join(["data", channel_name, period_name]) + "/"
        if file_path_postfix != "":
            path += file_path_postfix + "__/"
        path += file_name
        return path

    def get_mc_path(self, period_name, file_label, file_path_postfix=""):
        file_name = self.input_files["mc"][file_label]
        path = self.base_dir + "/".join(["mc", period_name]) + "/"
        if file_path_postfix != "":
            path += file_path_postfix + "__/"
        path += file_name
        return path

    def make_data_group(self, period_name, channel_name, *file_labels,
                        group_name="Data", file_path_postfix="", hist_prefix="",
                        hist_path_postfix=""):

        # use all files in a period
        if len(file_labels) == 0:
            file_labels = self.input_files["data"][channel_name][period_name].keys()

        # create FileGroup object
        # dictionary of files
        file_dict = dict()
        for file_label in file_labels:
            file_dict[file_label] = self.get_data_path(period_name, channel_name, file_label, file_path_postfix)

        sys_file_dict = None
        if self.sys_on:
            sys_file_dict = dict()
            for file_label in file_labels:
                sys_file_dict[file_label] = self.get_data_path(period_name, channel_name,
                                                               file_label, file_path_postfix="SYS")
        group = FileGroup(self.experiment_name, period_name, channel_name, group_name,
                          "data", file_dict, sys_file_dict, hist_path_postfix, hist_prefix)
        group.set_is_mc(False)
        return group

    def make_mc_group(self, period_name, channel_name, *file_labels,
                      group_name="mc", file_path_postfix="", hist_prefix="",
                      hist_path_postfix=""):

        # create FileGroup object
        # dictionary of files
        file_dict = dict()
        for file_label in file_labels:
            file_dict[file_label] = self.get_mc_path(period_name, file_label, file_path_postfix)

        sys_file_dict = None
        if self.sys_on:
            sys_file_dict = dict()
            for file_label in file_labels:
                sys_file_dict[file_label] = self.get_mc_path(period_name, file_label, file_path_postfix="SYS")

        group = FileGroup(self.experiment_name, period_name, channel_name, group_name,
                          "mc", file_dict, sys_file_dict, hist_path_postfix, hist_prefix)
        return group
