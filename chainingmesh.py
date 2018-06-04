

import numpy as np


def ChainingMesh3D:
    def __init__(self,cell_size,periodic=True):
        self.periodic=periodic
        self.cell_size = cell_size
    
    def set_data(self,x,y,z,min_val=None,max_val=None):
        if(min_val==None):
            min_val = min(min(x),min(y),min(z))
        if(max_val==None):
            max_val = max(max(x),max(y),max(z))
        self.cell_length = (max_val-min_val)/cell_size

    
    def get_index(self,i,j,k):
        
