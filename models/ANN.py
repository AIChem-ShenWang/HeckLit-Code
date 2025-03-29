import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from .fds import *

# 5-Layer ANN Regressor
class ANN(nn.Module):
    def __init__(self, input_size, FDS=None):
        super(ANN, self).__init__()
        self.input_size = input_size
        self.encoder = nn.Sequential(
            nn.Linear(input_size, 1000),
            nn.ReLU(),
            nn.Linear(1000, 5000),
            nn.ReLU(),
            nn.BatchNorm1d(5000),
        )
        self.decoder = nn.Sequential(
            nn.Linear(5000, 1000),
            nn.ReLU(),
            nn.BatchNorm1d(1000),
            nn.Linear(1000, 100),
            nn.ReLU(),
            nn.Linear(100, 1)
        )
        self.FDS = FDS

    def forward(self, x, y, epoch, device):
      x = self.encoder(x)
      if self.FDS is not None:
            if epoch >= 1:
                self.FDS.cpu()
                x = self.FDS.smooth(x.cpu(), y.cpu(), epoch).to(device)
                y = y.to(device)
      x = self.decoder(x.to(device))
      return x

# Model Evaluation Function
def RMSE(pred,true):
    diff_2 = (pred - true)**2
    return np.sqrt(diff_2.mean())

def R2(pred, true):
    u = ((true - pred) ** 2).sum()
    v = ((true - true.mean()) ** 2).sum()
    r2 = 1 - u / v
    return r2