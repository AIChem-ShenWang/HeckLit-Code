import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

# 5-Layer ANN Regressor
class ANN(nn.Module):
    def __init__(self, input_size):
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


    def forward(self, x):
      x = self.encoder(x)
      x = self.decoder(x)
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