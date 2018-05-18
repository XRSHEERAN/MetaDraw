class Common_Config(object):
    def __init__(self):
        self.time_step = 19 #at most time_step+1 length peptides
        self.max_ion_charge = 2
        self.ion_types = ['y','b']
        self.ion_terms = {'b':'n','y':'c','c':'n','z':'c','b_mod_loss':'n','y_mod_loss':'c'}
        self.from_CTerm = False
        self.modini = r'e:\GitHub\pdeep\config_file\modification.ini'
        self.fixmod = ['Carbamidomethyl[C]']
        self.varmod = []
        self.var_mod_num = 0
        self.mod_neutral_loss = False
        self.fragmentation = 'HCD'
        

class HCD_test_yion_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['y']
        
class HCD_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['y','b']

class HCD_ProteomeTools_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['y','b']
        self.fixmod = []

class ETD_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['c','z']
        self.fragmentation = 'ETD'

class ETD_ProteomeTools_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['c','z']
        self.fixmod = []
        self.fragmentation = 'ETD'
        
class EThcD_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['b','y','c','z']
        self.fragmentation = 'EThcD'
        
class EThcD_ProteomeTools_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['b','y','c','z']
        self.fixmod = []
        self.fragmentation = 'EThcD'

class HCD_pho_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['y','b']
        self.varmod = ['Phospho[S]','Phospho[T]','Phospho[Y]']
        self.var_mod_num = 1
        
class ETD_pho_Config(Common_Config):
    def __init__(self):
        super(self.__class__, self).__init__()
        self.ion_types = ['c','z']
        self.varmod = ['Phospho[S]','Phospho[T]','Phospho[Y]']
        self.var_mod_num = 1
        self.fragmentation = 'ETD'
        
