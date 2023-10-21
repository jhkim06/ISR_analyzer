from UnFolder import UnFolder
from ISRHists import ISRHists
import re


class ISRUnFolders:
    def __init__(self, channel, period, mass_edges,
                 file_path_postfix='', hist_path_postfix='', is_2d=False):

        self.channel = channel  # "ee2016a"
        self.period = period
        self.mass_edges = mass_edges
        self.file_path_postfix = file_path_postfix
        self.hist_path_postfix = hist_path_postfix
        self.is_2d = is_2d

        self.mass_unfolders = []
        self.pt_unfolders = []

        self.folded_isr_hists = None
        self.unfolded_isr_hists = None
        self.full_phase_isr_hists = None

    def set_unfolders(self,  input_hist, matrix, bg_hist=None, unfolder_name='mass'):
        if unfolder_name == "mass":
            unfolder_container = self.mass_unfolders
        else:
            unfolder_container = self.pt_unfolders

        if self.is_2d:
            # UnFolder will extract 1D from 2D input
            unfolder_container.append(UnFolder(matrix, input_hist, bg_hist))
        else:
            for i in range(len(input_hist)):
                if bg_hist is None:
                    unfolder_container.append(UnFolder(matrix[i], input_hist[i]))
                else:
                    unfolder_container.append(UnFolder(matrix[i], input_hist[i], bg_hist[i]))

    def unfold(self):
        for i in range(len(self.mass_unfolders)):
            self.mass_unfolders[i].do_unfold()

        folded_pt = []
        unfolded_pt = []
        for i in range(len(self.pt_unfolders)):
            self.pt_unfolders[i].do_unfold()

            folded_pt.append(self.pt_unfolders[i].get_folded_hist())
            unfolded_pt.append(self.pt_unfolders[i].get_unfolded_hist())

        unfolded_mass = self.mass_unfolders[0].get_unfolded_hist()
        self.unfolded_isr_hists = ISRHists(unfolded_mass, *unfolded_pt,
                                           mass_windows=self.mass_edges)
        folded_mass = self.mass_unfolders[0].get_folded_hist()
        self.folded_isr_hists = ISRHists(folded_mass, *folded_pt,
                                         mass_windows=self.mass_edges)

    def acceptance_correction(self, hist_producer, unfolder_name='mass'):
        if unfolder_name == 'mass':
            unfolder = self.mass_unfolders
        else:
            unfolder = self.pt_unfolders

        full_phase_data_hist = []
        for i in range(len(unfolder)):
            unfolded_data_hist = unfolder[i].get_unfolded_hist()

            full_phase_hist_name = re.sub(r'(__)', '_acceptance__',
                                          unfolded_data_hist.hist_name)
            unfolded_data_hist.set_hist_name(full_phase_hist_name)

            unfolded_mc_hist = unfolder[i].get_expectation_hist(use_matrix=True)
            full_phase_hist = hist_producer.get_signal_hist(full_phase_hist_name)
            acceptance = full_phase_hist.divide(unfolded_mc_hist)

            full_phase_data_hist.append(unfolded_data_hist.multiply(acceptance))
        return full_phase_data_hist

    def acceptance_corrections(self, hist_producer):
        hist_producer.set_base_configs(self.period, self.channel, self.file_path_postfix, self.hist_path_postfix)
        full_phase_mass_data_hist = self.acceptance_correction(hist_producer)
        full_phase_pt_data_hist = self.acceptance_correction(hist_producer, 'pt')

        self.full_phase_isr_hists = ISRHists(full_phase_mass_data_hist[0], *full_phase_pt_data_hist,
                                             mass_windows=self.mass_edges)

    def make_mc_isr_hists(self, unfolded=False, use_matrix=True):
        print("make mc isr hists")
        mass_mc_hist = self.mass_unfolders[0].get_expectation_hist(unfolded=unfolded, use_matrix=use_matrix)
        pt_mc_hists = []
        for i in range(len(self.pt_unfolders)):
            pt_mc_hists.append(self.pt_unfolders[i].get_expectation_hist(unfolded=unfolded, use_matrix=use_matrix))
        out_isr_hists = ISRHists(mass_mc_hist, *pt_mc_hists,
                                 mass_windows=self.mass_edges)
        return out_isr_hists

    def get_mc_full_phase_isr_hists(self, hist_producer):
        hist_producer.set_base_configs(self.period, self.channel,
                                       hist_path_postfix=self.hist_path_postfix)

        mass_hist = self.full_phase_isr_hists.get_mass_hist()

        mass_hist_name = mass_hist.hist_name
        if self.is_2d:
            matches = re.findall(r'(\[[^]]*])', mass_hist_name)
            mass_hist_name = matches[0] + "_" + matches[1] + "_" + matches[2]
        hist_producer.set_configs(mass_hist_name)

        mc_mass_hist = hist_producer.get_total_expectation_hist(exp_type='signal')

        pt_hists = self.full_phase_isr_hists.get_pt_hist(-1)
        mc_pt_hists = []
        for i in range(len(pt_hists)):
            pt_hist_name = pt_hists[i].hist_name
            if self.is_2d:
                matches = re.findall(r'(\[[^]]*])', pt_hist_name)
                pt_hist_name = matches[0] + "_" + matches[1] + "_" + matches[2]
            hist_producer.set_configs(pt_hist_name)

            mc_pt_hists.append(hist_producer.get_total_expectation_hist(exp_type='signal'))

        out_isr_hists = ISRHists(mc_mass_hist, *mc_pt_hists,
                                 mass_windows=self.mass_edges)
        return out_isr_hists

    def get_isr_hists(self, level_name='folded'):
        if level_name == 'folded':
            return self.get_folded_isr_hists()
        elif level_name == 'unfolded':
            return self.get_unfolded_isr_hists()
        elif level_name == 'full_phase':
            return self.get_full_phase_isr_hists()
        else:
            print('check ' + level_name)
            return

    def get_mc_isr_hists(self, level_name='folded', hist_producer=None):
        if level_name == 'folded':
            return self.make_mc_isr_hists(unfolded=False)
        elif level_name == 'unfolded':
            return self.make_mc_isr_hists(unfolded=True)
        elif level_name == 'full_phase':
            return self.get_mc_full_phase_isr_hists(hist_producer)
        else:
            print('check ' + level_name)
            return

    def get_folded_isr_hists(self):
        return self.folded_isr_hists

    def get_unfolded_isr_hists(self):
        return self.unfolded_isr_hists

    def get_full_phase_isr_hists(self):
        return self.full_phase_isr_hists
