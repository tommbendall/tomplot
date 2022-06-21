"""
Declares a class for holding information about how to compute a diagnostic field
"""

class DiagnosticInfo(object):

    def __init__(self, name, prognostic_variables, evaluator):

        self.name = name
        self.is_data_extracted = False
        self.prognostic_variables = prognostic_variables
        self.evaluator = evaluator


    def extract_lfric_data(self, data_file, time_idx):

        self.data = {}

        for prognostic in self.prognostic_variables:
            self.data[prognostic] = data_file[prognostic][time_idx,:,:]

        self.is_data_extracted = True

    def extract_data(self, data_file, time_idx):

        raise NotImplementedError('extract_data is not implemented')

    
    def evaluate(self):

        if not self.is_data_extracted:
            raise ValueError('Before diagnostic can be evaluated, data must be extracted')
        
        return self.evaluator(self.data)
