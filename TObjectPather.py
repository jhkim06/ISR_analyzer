

class TObjectPather:
    def __init__(self, period_name, channel_name, hist_path_postfix, hist_prefix):
        self.channel_name = channel_name
        self.period_name = period_name
        self.hist_path_postfix = hist_path_postfix

        if hist_prefix != "":
            self.hist_prefix = hist_prefix + "_"
        else:
            self.hist_prefix = ""

    def get_hist_path(self, sys=False):
        return self.get_path(sys) + self.hist_prefix

    def get_path(self, sys=False):
        if sys or self.hist_path_postfix == "":
            return self.channel_name + self.period_name + "/"
        else:
            return self.channel_name + self.period_name + "/" + self.hist_path_postfix + "/"
