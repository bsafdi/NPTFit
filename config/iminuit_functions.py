import numpy

class fake_code:
    def __init__(self, n,var_names=['False']):
        self.co_argcount = n
        if var_names[0] != 'False':
            self.co_varnames = tuple(var_names)
        else:
            self.co_varnames = tuple(map(str, range(n)))
 
 ######
 ###input: n is number of model parameters, ll is the likelihood function that takes as input 
 #####ll([x1,x2,...,xn])
 ####convert ll function into form necessary for iminuit
class call_ll:
    def __init__(self, n,ll,var_names=['False']):
        self.n = n
        self.func_code = fake_code(n,var_names=var_names)
        self.ll = ll
    def __call__(self, *args):
        return - self.ll([arg for arg in args])