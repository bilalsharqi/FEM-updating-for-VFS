import numpy as np
import os
import platform
import shutil
from openmdao.components.external_code_comp import ExternalCodeComp

nastran = '/apps/global/msc/nastran/2012.2/rhel5/x86_64/bin/nast20122'
scratch_dir = '/home/deatonjd/scratch'

this_file_dir = os.path.abspath(os.path.dirname(__file__))
test_dir = os.path.abspath(os.path.join(this_file_dir, 'test'))

n_input = 287  # 287 panel thicknesses
n_output = 56*6*6  # (56 beam nodes)(6 degrees-of-freedom)(6 load cases)=2016


class NastranSensExternalCodeComponent(ExternalCodeComp):

    # def initialize(self):
    #     # Add in parameterized loads if necessary
    #     self.options.declare('loads', types=[])

    def setup(self):
        self.add_input('thickness', 0.010*np.ones((n_input,)), units='m')
        self.add_output('disp', np.zeros((n_output,)), units='m')

        # self.sol101_f06 = os.path.abspath('sol101.f06')
        # self.sol200_sens_csv = os.path.abspath('sol200_sens.csv')

        # Below is a dummy command. Command is changed between function
        # evaluation [.compute()] or gradient evaluation [.compute_partials()]
        # since we need to run separate Nastran input decks.
        self.options['command'] = ['pwd']

        self.declare_partials(of='*', wrt='*')

    def compute(self, inputs, outputs):

        # Cleanup files from any old runs.
        self._cleanup_files([os.path.abspath('uCRM-135_mat_and_props.bdf'),
                             os.path.abspath('sol101.f06'),
                             os.path.abspath('sol101.DBALL'),
                             os.path.abspath('sol101.MASTER'),
                             os.path.abspath('sol101.f04'),
                             os.path.abspath('sol101.IFPDAT'),
                             os.path.abspath('sol101.log'),
                             os.path.abspath('sol101.op2'),
                             os.path.abspath('sol101.xdb')])

        # Produce input file.
        thickness = inputs['thickness']
        input_file_data = dict()
        for i in range(0, n_input):
            input_file_data['xE{}x'.format(i+1)] = 72.0e9
            input_file_data['xNU{}x'.format(i+1)] = 0.33
            input_file_data['xRHO{}x'.format(i+1)] = 2780.0
            input_file_data['xT{}x'.format(i+1)] = thickness[i]
        self._replace_keywords(input_file_data,
                               os.path.abspath('uCRM-135_mat_and_props.template'),
                               os.path.abspath('uCRM-135_mat_and_props.bdf'))

        # Set required input/output files.
        self.options['external_input_files'] = [os.path.abspath('sol101.bdf'),
                                                os.path.abspath('uCRM-135_beam_forces.bdf'),
                                                os.path.abspath('uCRM-135_beam_full_rbe3_equal_weight.bdf'),
                                                os.path.abspath('uCRM-135_beam_grids.bdf'),
                                                os.path.abspath('uCRM-135_mat_and_props.bdf'),
                                                os.path.abspath('uCRM-135_spc.bdf'),
                                                os.path.abspath('uCRM-135_wingbox_coarse_mesh.bdf')]
        self.options['external_output_files'] = [os.path.abspath('sol101.f06')]

        # Run Nastran SOL101 input deck to calculate displacements.
        self.options['command'] = [nastran, 'sol101.bdf', 'batch=no',
                                   'sdir='+scratch_dir, 'old=no', 'news=no']
        if 'Darwin' in platform.platform():
            # Don't try to run Nastran on Mac, we're testing if on Mac.
            shutil.copyfile(os.path.join(test_dir, 'sol101.f06'),
                            os.path.join(this_file_dir, 'sol101.f06'))
        else:
            super(NastranSensExternalCodeComponent, self).compute(inputs, outputs)

        # Parse/read output files and populate outputs. Must unwrap list & dict
        # from packaged output from parser methods.
        disp_list = []
        for subcase in self._get_f06_disp(os.path.abspath('sol101.f06')):
            for node_id, disps in subcase.items():
                for d in disps:
                    disp_list.append(d)
        outputs['disp'] = np.array(disp_list)
        # outputs['disp'] = 0.111*np.ones((n_output,))

    def compute_partials(self, inputs, partials):

        # Cleanup files from any old runs.
        self._cleanup_files([os.path.abspath('uCRM-135_design_data.bdf'),
                             os.path.abspath('uCRM-135_mat_and_props.bdf'),
                             os.path.abspath('sol200_sens.csv'),
                             os.path.abspath('sol200.DBALL'),
                             os.path.abspath('sol200.f04'),
                             os.path.abspath('sol200.f06'),
                             os.path.abspath('sol200.IFPDAT'),
                             os.path.abspath('sol200.log'),
                             os.path.abspath('sol200.MASTER'),
                             os.path.abspath('sol200.op2'),
                             os.path.abspath('sol200.xdb')])

        # Produce input files.
        thickness = inputs['thickness']
        input_file_data = dict()
        for i in range(0, n_input):
            input_file_data['xE{}x'.format(i + 1)] = 72.0e9
            input_file_data['xNU{}x'.format(i + 1)] = 0.33
            input_file_data['xRHO{}x'.format(i + 1)] = 2780.0
            input_file_data['xT{}x'.format(i + 1)] = thickness[i]
        input_file_data['xDV_LBx'] = 1.e-6
        input_file_data['xDV_UBx'] = 1.0
        self._replace_keywords(input_file_data,
                               os.path.abspath('uCRM-135_mat_and_props.template'),
                               os.path.abspath('uCRM-135_mat_and_props.bdf'))
        self._replace_keywords(input_file_data,
                               os.path.abspath('uCRM-135_design_data.template'),
                               os.path.abspath('uCRM-135_design_data.bdf'))

        # Set required input/output files.
        self.options['external_input_files'] = [os.path.abspath('sol200.bdf'),
                                                os.path.abspath('uCRM-135_beam_forces.bdf'),
                                                os.path.abspath('uCRM-135_beam_full_rbe3_equal_weight.bdf'),
                                                os.path.abspath('uCRM-135_beam_grids.bdf'),
                                                os.path.abspath('uCRM-135_design_data.bdf'),
                                                os.path.abspath('uCRM-135_mat_and_props.bdf'),
                                                os.path.abspath('uCRM-135_spc.bdf'),
                                                os.path.abspath('uCRM-135_wingbox_coarse_mesh.bdf')]
        self.options['external_output_files'] = [os.path.abspath('sol200.f06'),
                                                 os.path.abspath('sol200_sens.csv')]

        # Run Nastran SOL200 input deck to calculate displacement sensitivities.
        self.options['command'] = [nastran, 'sol200.bdf', 'batch=no',
                                   'sdir=' + scratch_dir, 'old=no', 'news=no']
        outputs = {}
        if 'Darwin' in platform.platform():
            # Don't try to run Nastran on Mac, we're testing if on Mac.
            shutil.copyfile(os.path.join(test_dir, 'sol200.f06'),
                            os.path.join(this_file_dir, 'sol200.f06'))
            shutil.copyfile(os.path.join(test_dir, 'sol200_sens.csv'),
                            os.path.join(this_file_dir, 'sol200_sens.csv'))
        else:
            super(NastranSensExternalCodeComponent, self).compute(inputs, outputs)

        # Parse output and setup partials.
        du_dx = self._get_sens_disp(os.path.abspath('sol200_sens.csv'))
        partials['disp', 'thickness'] = self._organize_partials(du_dx)

    #
    # Methods below here are helpers for writing input files, parsing output
    # files, etc. We make them class/static methods since we will not change
    # object state in any of them.
    #

    @classmethod
    def _cleanup_files(cls, files_to_delete):
        """Delete specified files if they exist."""
        for file in files_to_delete:
            if os.path.exists(file):
                os.remove(file)

    @classmethod
    def _replace_keywords(cls, keywords, template_file, new_file):
        """Replaces keywords in a template to produce a new file.

        Parameters:
            keywords (dict): Dictionary of keyword, value pairs.
            template_file (str): Path to a file templated with keywords.
            new_file: Path to file to be created.

        Returns:
            None
        """
        with open(template_file, 'r') as template_f:
            with open(new_file, 'w') as new_f:
                for line in template_f:
                    for (key, value) in keywords.items():
                        # Not the most ideal way to do this but gets the job
                        # done for now.
                        if 'E' in key:
                            line = line.replace(key, '{0:.2e}'.format(value))
                        elif 'NU' in key:
                            line = line.replace(key, '{0:.2f}'.format(value))
                        elif 'RHO' in key:
                            line = line.replace(key, '{0:.1f}'.format(value))
                        elif 'T' in key:
                            line = line.replace(key, '{0:.6f}'.format(value))
                        elif 'DV_LB' in key:
                            line = line.replace(key, '{0:.6f}'.format(value))
                        elif 'DV_UB' in key:
                            line = line.replace(key, '{0:.6f}'.format(value))
                        else:
                            line = line.replace(key, 'XX')
                    new_f.write(line)

    @classmethod
    def _get_f06_disp(cls, f06_file):
        """Retrieves displacement from Nastran .f06 output file.

        Parses an .f06 to read the displacement for multiple SOL101 subcases.

        Parameters:
            f06_file (str): path to a Nastran .f06 output file.
        Returns:
            list(dict): list (subcase) of dict (node_id) containing displacement
                        vector for each node. Each top level list element is a
                        subcase in the order found in .f06 file.
        """

        disp = []
        with open(f06_file, "r") as f:
            for line in f:
                if "D I S P L A C E M E N T   V E C T O R" in line:
                    disp.append(cls._parse_displacements(f))
        return disp

    @staticmethod
    def _parse_displacements(f):
        """Helper function to extract formatted displacement data.

        We only extract data corresponding to specified nodes after the caller
        has located a section of displacement data in the file object.
        """
        disp = dict()
        nodes = range(30000001, 30000057)
        for node in nodes:
            line = f.readline()
            while str(node) not in line:
                line = f.readline()
            disp[node] = [float(val) for val in line.split()[2:]]
        return disp

    @classmethod
    def _get_sens_disp(cls, sens_file):
        """Retrieves sensitivities from Nastran SOL200 .csv output file.

        Parses output file to get sensitivity data corresponding to nodal
        displacement (multiple components) for multiple load subcases.

        Parameters:
            sens_file (str): path to a Nastran .csv output file from SOL200.
        Returns:
            dict: nested dictionary organized by subcase id [1-6], node
                  id [30000001-30000056], and component id [1-6] as
                  du_dx[subcase_id][node_id][comp_id] = [list of du_dxi].
        """
        # Top level index = load case (subcase id).
        du_dx = {1: dict(), 2: dict(), 3: dict(),
                 4: dict(), 5: dict(), 6: dict()}
        with open(sens_file, "r") as f:
            for line in f:
                if 'Design Response ID' in line:
                    line1 = f.readline().split(',')  # Response/node/component/case identifier line.
                    f.readline()                     # DV_ID line (not used).
                    line3 = f.readline().split(',')  # DV value line.
                    if line1[3] == "            ":
                        # this condition skips the first weight response in
                        # the csv output file.
                        pass
                    else:
                        grid_id = int(line1[3])
                        comp_id = int(line1[4])
                        subcase = int(line1[6])
                        if grid_id not in du_dx[subcase]:
                            du_dx[subcase][grid_id] = dict()
                        du_dx[subcase][grid_id][comp_id] = [float(du_dxi) for du_dxi in line3]
        return du_dx

    @classmethod
    def _organize_partials(cls, sens_dict):
        """Helper function to reshape sensitivity output.

        This function reshapes the organized sensitivity output in a nested
        dictionary organized as:
               du_dx[subcase_id][node_id][comp_id] = [list of du_dxi]
        into a numpy array organized as:
               |du(load1_grid1_comp1)_dx1, du(load1_grid1_comp1)_dx2, ...
               |du(load1_grid1_comp2)_dx1, du(load1_grid1_comp2)_dx2, ...
               |          :                          :
               |du(load1_grid2_comp1)_dx1, du(load1_grid2_comp1)_dx2, ...
               |du(load1_grid2_comp2)_dx1, du(load1_grid2_comp2)_dx2, ...
               |          :                          :
               |du(load2_grid1_comp1)_dx1, du(load1_grid1_comp1)_dx2, ...
               |du(load2_grid1_comp2)_dx1, du(load1_grid1_comp2)_dx2, ...
               |          :                          :
        Parameters:
            sens_dict (dict): Nested dictionary of sensitivity data.
        Returns:
            np.array: (n_disp_responses)x(n_dvs) size array of partials.
        """
        partials = []
        # Use sorted() to safely ensure we get things in order regardless of
        # version of Python (rather than changing to OrderedDict or expecting
        # Python 3)
        for s_id in sorted(sens_dict.keys()):
            for g_id in sorted(sens_dict[s_id].keys()):
                for c_id in sorted(sens_dict[s_id][g_id].keys()):
                    partials.append(sens_dict[s_id][g_id][c_id])
        return partials


def reference_thicknesses():
    """Reference thickness from uCRM-135 as downloaded from MDOLAB website. """
    return np.array([0.0131, 0.0118, 0.00946, 0.00785, 0.0111, 0.00726, 0.00765,
                     0.00853, 0.00816, 0.00752, 0.00711, 0.00677, 0.00603,
                     0.00571, 0.00615, 0.00728, 0.0066, 0.00632, 0.00894,
                     0.00615, 0.00594, 0.00583, 0.00581, 0.00584, 0.00593,
                     0.00603, 0.00608, 0.00613, 0.0061, 0.00603, 0.00593,
                     0.00588, 0.00571, 0.00589, 0.00728, 0.0111, 0.00926,
                     0.0109, 0.0121, 0.0111, 0.00728, 0.00546, 0.00528, 0.00528,
                     0.00528, 0.00528, 0.00528, 0.00528, 0.00528, 0.00528,
                     0.00528, 0.00528, 0.00528, 0.00528, 0.00528, 0.00528,
                     0.00703, 0.00637, 0.00662, 0.00943, 0.0133, 0.0146, 0.0108,
                     0.00711, 0.00728, 0.00528, 0.00528, 0.00542, 0.00579,
                     0.00616, 0.00698, 0.00847, 0.0123, 0.0118, 0.00793,
                     0.00725, 0.00584, 0.00728, 0.00652, 0.00633, 0.00929,
                     0.00728, 0.0111, 0.00728, 0.00879, 0.00607, 0.00728,
                     0.00584, 0.00575, 0.00562, 0.00565, 0.00569, 0.00749,
                     0.0056, 0.00552, 0.00728, 0.00538, 0.0053, 0.00528,
                     0.00528, 0.00528, 0.00528, 0.00528, 0.00528, 0.00528,
                     0.00528, 0.00528, 0.00528, 0.00528, 0.00528, 0.0057,
                     0.00502, 0.00572, 0.00607, 0.00487, 0.00487, 0.00487,
                     0.00496, 0.00511, 0.00513, 0.00511, 0.00516, 0.00509,
                     0.00492, 0.00487, 0.00487, 0.00664, 0.00531, 0.00507,
                     0.00493, 0.00487, 0.00487, 0.00487, 0.00487, 0.00487,
                     0.00487, 0.00487, 0.00487, 0.00487, 0.00487, 0.00487,
                     0.00487, 0.00487, 0.00487, 0.00487, 0.00487, 0.00487,
                     0.00487, 0.00487, 0.00487, 0.00487, 0.00487, 0.00487,
                     0.00487, 0.00487, 0.00487, 0.00487, 0.00487, 0.00487,
                     0.00487, 0.00487, 0.00487, 0.00487, 0.00487, 0.00487,
                     0.00487, 0.00487, 0.0046, 0.00466, 0.00466, 0.00464,
                     0.00518, 0.00463, 0.00569, 0.00461, 0.0062, 0.00498,
                     0.00675, 0.00546, 0.00745, 0.00594, 0.00813, 0.00643,
                     0.00876, 0.00717, 0.00938, 0.00795, 0.00998, 0.00873,
                     0.0106, 0.00949, 0.0112, 0.0103, 0.0119, 0.011, 0.0127,
                     0.0118, 0.0136, 0.0125, 0.0145, 0.0132, 0.0154, 0.014,
                     0.0163, 0.0149, 0.0172, 0.0157, 0.0181, 0.0165, 0.0189,
                     0.0174, 0.0198, 0.0182, 0.0206, 0.0191, 0.0215, 0.02,
                     0.0223, 0.0209, 0.0232, 0.0217, 0.0241, 0.0226, 0.025,
                     0.0235, 0.0259, 0.0244, 0.0264, 0.0252, 0.0271, 0.0259,
                     0.028, 0.0268, 0.0288, 0.0275, 0.0286, 0.0274, 0.0285,
                     0.028, 0.0256, 0.035, 0.0262, 0.0352, 0.0261, 0.0342,
                     0.0255, 0.0333, 0.0247, 0.0322, 0.0239, 0.0311, 0.0231,
                     0.03, 0.0223, 0.0289, 0.0215, 0.0278, 0.0207, 0.0267,
                     0.0199, 0.0256, 0.0191, 0.0245, 0.0184, 0.0245, 0.0177,
                     0.0234, 0.0179, 0.0223, 0.0171, 0.0213, 0.0173, 0.0268,
                     0.0185, 0.0257, 0.0178, 0.0246, 0.017, 0.0235, 0.0164,
                     0.0339, 0.0115, 0.00487, 0.0164, 0.045, 0.0141, 0.0142])


# Test component.
if __name__ == '__main__':
    from openmdao.api import Problem, IndepVarComp

    prob = Problem()
    model = prob.model

    # Initialize independent variable component with reference thicknesses.
    model.add_subsystem('des_vars', IndepVarComp('x', reference_thicknesses()))

    # Keep reference to nastran components so we can test class methods.
    nastran_component = NastranSensExternalCodeComponent()
    model.add_subsystem('nastran', nastran_component)

    model.connect('des_vars.x', 'nastran.thickness')

    prob.setup()

    prob.run_model()
    print(prob['nastran.disp'])

    # Force execution of .compute() with initialized inputs
    #nastran_component.compute(nastran_component._inputs, {})

    # Force execution of .compute_partials() with initialized inputs
    nastran_component.compute_partials(nastran_component._inputs, {})
    du_dx = nastran_component._get_sens_disp("sol200_sens.csv")
    du_dx_flat = nastran_component._organize_partials(du_dx)

    # Run check partials.
    # data = prob.check_partials()



