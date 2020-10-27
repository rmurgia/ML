import numpy as np
import sys,os
import scipy.special as si

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils import data


# num_redshifts = 39
# num_los = 5000
# num_redshifts = 1
# num_los = 10
num_pixels = 2048
folder = "../ML_catalogues/NCDM_1/z=2.2"

# print(np.shape(x))

# data loader
class Dataset(data.Dataset):
    def __init__(self):
        rhoker_H, rhoker_H1 tempker_H1, velker_H1 = np.loadtxt()

    #     y = np.zeros(300, dtype=np.float32)
    #     indexes = np.where(x<1.0)[0]
    #     y[indexes] = np.sin((x[indexes]-1.0)*30.0)
    #     indexes = np.where((x>=1.0) & (x<1.5))[0]
    #     y[indexes] = (x[indexes])**2-1.0
    #     indexes = np.where((x>=1.5))[0]
    #     y[indexes] = si.jv(0,(x[indexes]-1.5)*40)*2.0 -0.75

    #     self.x = torch.from_numpy(x)
    #     self.y = torch.from_numpy(y)
    #     self.len = self.x.shape[0]

    # def __len__(self):
    #     return self.len

    # def __getitem__(self,index):
    #     return self.x[index], torch.tensor([self.x[index], self.x[index]**2]), self.y[index]