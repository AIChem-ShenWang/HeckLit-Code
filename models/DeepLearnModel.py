import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from .fds import *
import torch_geometric.nn as pyg_nn
from torch_geometric.data import Data
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

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

    def forward(self, x, y=None, epoch=None, device=None):
      x = self.encoder(x)
      if self.FDS is not None:
            if epoch >= 1:
                self.FDS.cpu()
                x = self.FDS.smooth(x.cpu(), y.cpu(), epoch).to(device)
                y = y.to(device)
      x = self.decoder(x.to(device))
      return x

# MMHRP
class GATU(nn.Module):
    def __init__(self, node_feature_num, channels, heads):
        super(GATU, self).__init__()
        self.conv1 = pyg_nn.GATConv(node_feature_num, channels[0], heads=heads)
        self.norm1 = nn.BatchNorm1d(channels[0] * heads)
        self.conv2 = pyg_nn.GATConv(channels[0] * heads, channels[1], heads=heads)

    def forward(self, data):
        x = data["x"].clone().detach().float()
        edge_index = data["edge_index"].clone().detach()
        batch = data["batch"].clone().detach()

        # GAT
        x1 = self.conv1(x, edge_index)
        x1 = self.norm1(x1)
        x1 = F.relu(x1)

        x2 = self.conv2(x1, edge_index)
        x2 = F.relu(x2)
        x = torch.cat([x1, x2], 1)

        # Pooling
        x_mean = pyg_nn.global_mean_pool(x, batch=batch)
        x_max = pyg_nn.global_max_pool(x, batch=batch)
        x = torch.cat([x_mean, x_max], 1)

        return x

class GCNU(nn.Module):
    def __init__(self, node_feature_num, channels):
        super(GCNU, self).__init__()
        self.conv1 = pyg_nn.GCNConv(node_feature_num, channels[0])
        self.norm1 = nn.BatchNorm1d(channels[0])
        self.conv2 = pyg_nn.GCNConv(channels[0], channels[1])

    def forward(self, data):
        x = data["x"].clone().detach().float()
        edge_index = data["edge_index"].clone().detach()
        batch = data["batch"].clone().detach()

        # GAT
        x1 = self.conv1(x, edge_index)
        x1 = self.norm1(x1)
        x1 = F.relu(x1)

        x2 = self.conv2(x1, edge_index)
        x2 = F.relu(x2)
        x = torch.cat([x1, x2], 1)

        # Pooling
        x_mean = pyg_nn.global_mean_pool(x, batch=batch)
        x_max = pyg_nn.global_max_pool(x, batch=batch)
        x = torch.cat([x_mean, x_max], 1)

        return x

class BiGRU(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, output_size, device):
        super(BiGRU, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.gru = nn.GRU(
            input_size=input_size,
            hidden_size=hidden_size,
            num_layers=num_layers,
            batch_first=True,  # 输入时seq_len与batch_size顺序颠倒
            bidirectional=True
        )
        self.fc = nn.Linear(2 * hidden_size, output_size)
        self.device = device

    def init_hidden(self, x):
      return torch.zeros(self.num_layers * 2, x.size(0), self.hidden_size).to(self.device)

    def forward(self, x):
        h0 = self.init_hidden(x)
        out, _ = self.gru(x, h0) # x is input, size (batch, seq_len, input_size)
        out = out[:, -1, :]
        out = self.fc(out)
        return out

class MMHRP_GCL(nn.Module):
    def __init__(self,
                 GraphEncoder: dict,
                 TextEncoder: dict,
                 ModalityAlignment: dict,
                 Decoder: dict,
                 emb_size: int = 128, # embedding vector size for each encoder
                 device: str = torch.device('cuda'),
                 FDS=None
                 ):
        super(MMHRP_GCL, self).__init__()

        self.device = device
        self.GraphEncoder = GraphEncoder
        self.TextEncoder = TextEncoder
        self.FDS = FDS

        # 1.Encoder Parser
        # 1.1 GraphEncoder Parser
        Graph_params = ["NodeFeatNum", "Channels", "Heads"]
        for key in GraphEncoder.keys():
            if key not in Graph_params:
                raise Exception("%s is not the param in Graph Encoder")

            if key == "NodeFeatNum":
                NodeFeatNum = GraphEncoder[key]  # 8
            if key == "Channels":
                GAT_Channels = GraphEncoder[key]  # [32, 64]
                assert len(GAT_Channels) == 2
            if key == "Heads":
                GAT_Heads = GraphEncoder[key]  # 4

        # GATU output size
        GATU_OutSize = 0
        for i in GAT_Channels:
            GATU_OutSize += i * GAT_Heads
        self.GATU_OutSize = GATU_OutSize * 2

        # GraphEncoder
        self.GATU_ReaPro = GATU(NodeFeatNum, GAT_Channels, GAT_Heads)
        self.GATU_CatSol = GATU(NodeFeatNum, GAT_Channels, GAT_Heads)
        # GATU for Reactants and Products
        self.ReaProEncoder = nn.Sequential(
            self.GATU_ReaPro,
            nn.Linear(self.GATU_OutSize, emb_size),
        )
        # GATU for Catalysts and Solvents
        self.CatSolEncoder = nn.Sequential(
            self.GATU_CatSol,
            nn.Linear(self.GATU_OutSize, emb_size),
        )

        # 1.2 TextEncoder Parser
        Text_params = ["SmiFeatNum", "Heads", "BigruChannels", "BigruNumlayer"]
        for key in TextEncoder.keys():
            if key not in Text_params:
                raise Exception("%s is not the param in Text Encoder")

            if key == "SmiFeatNum":
                self.SmiFeatNum = TextEncoder[key]  # 128
            if key == "Heads":
                Trans_Heads = TextEncoder[key]  # 4
            if key == "BigruChannels":
                BigruChannels = TextEncoder[key]  # [128, 128]
                assert len(BigruChannels) == 2
                bigru_input, bigru_hidden = BigruChannels
            if key == "BigruNumlayer":
                bigru_num_layers = TextEncoder[key]  # 2

        # TextEncoder
        self.trans = nn.TransformerEncoderLayer(d_model=self.SmiFeatNum, nhead=Trans_Heads, batch_first=True,
                                                norm_first=True)
        self.bigru = BiGRU(bigru_input, bigru_hidden, bigru_num_layers, emb_size, device)
        self.RxnSmiEncoder = nn.Sequential(
            self.trans,
            self.bigru,
        )

        # 3. Modality Alignment Parser
        MA_params = ["Heads"]
        for key in ModalityAlignment.keys():
            if key not in MA_params:
                raise Exception("%s is not the param in ModalityAlignment")

            if key == "Heads":
                MA_Heads = ModalityAlignment[key]  # 4

        self.MA = nn.TransformerEncoderLayer(d_model=emb_size * 3, nhead=MA_Heads)

        # 4. Decoder Parser
        Decoder_params = ["Channels"]
        for key in Decoder.keys():
            if key not in Decoder_params:
                raise Exception("%s is not the param in Decoder")

            if key == "Channels":
                Decoder_Channels = Decoder[key]  # [1000, 500, 100]
                assert len(Decoder_Channels) == 3

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(emb_size * 3, Decoder_Channels[0]),
            nn.Dropout(0.8),
            nn.ReLU(),
            nn.BatchNorm1d(Decoder_Channels[0]),
            

            nn.Linear(Decoder_Channels[0], Decoder_Channels[1]),
            nn.ReLU(),

            nn.Linear(Decoder_Channels[1], Decoder_Channels[2]),
            nn.ReLU(),

            nn.Linear(Decoder_Channels[2], 1)
        )

    def forward(self, x, y=None, epoch=None, device=None):
        # import data
        ReaPro_data, CatSol_data, RxnSmi = x
        # Graph Modality
        graph_emb = torch.cat([self.ReaProEncoder(ReaPro_data),
                                   self.CatSolEncoder(CatSol_data)], dim=1)
        # Text Modality
        text_embed = self.RxnSmiEncoder(RxnSmi)
        x = torch.cat([graph_emb, text_embed], dim=1)


        # Modality Alignment
        x = self.MA(x)

        # FDS operation
        if self.FDS is not None:
            if epoch >= 1:
                self.FDS.cpu()
                x = self.FDS.smooth(x.cpu(), y.cpu(), epoch).to(device)
                y = y.to(device)

        # Decoder
        x = self.decoder(x)

        return x

class MMHRP_GFP(nn.Module):
    def __init__(self,
                 GraphEncoder: dict,
                 FPEncoder: dict,
                 ModalityAlignment: dict,
                 Decoder: dict,
                 emb_size: int = 128, # embedding vector size for each encoder
                 device: str = torch.device('cuda'),
                 FDS=None
                 ):
        super(MMHRP_GFP, self).__init__()

        # device
        self.device = device
        self.GraphEncoder = GraphEncoder
        self.FPEncoder = FPEncoder
        self.FDS = FDS

        # 1.Encoder Parser
        # 1.1 GraphEncoder Parser
        Graph_params = ["Type", "NodeFeatNum", "Channels", "Heads"]
        for key in GraphEncoder.keys():
            if key not in Graph_params:
                raise Exception("%s is not the param in Graph Encoder")

            if key == "Type":
                t = GraphEncoder[key]
                type = ["GAT", "GCN"]
                if t not in type:
                    raise Exception("%s is not a Graph Encoder Type")
            if key == "NodeFeatNum":
                NodeFeatNum = GraphEncoder[key]  # 8
            if key == "Channels":
                GE_Channels = GraphEncoder[key]  # [32, 64]
                assert len(GE_Channels) == 2
            if key == "Heads":
                GAT_Heads = GraphEncoder[key]  # 4

        # GATU output size
        GE_OutSize = 0
        if t == "GAT":
            for i in GE_Channels:
                GE_OutSize += i * GAT_Heads
        self.GE_OutSize = GE_OutSize * 2

        # GraphEncoder
        if t == "GAT":
            self.ReaPro = GATU(NodeFeatNum, GE_Channels, GAT_Heads)
            self.CatSol = GATU(NodeFeatNum, GE_Channels, GAT_Heads)
        if t == "GCN":
            self.ReaPro = GCNU(NodeFeatNum, GE_Channels)
            self.CatSol = GCNU(NodeFeatNum, GE_Channels)
        # GATU for Reactants and Products
        self.ReaProEncoder = nn.Sequential(
            self.ReaPro,
            nn.Linear(self.GE_OutSize, int(emb_size/2))
        )
        # GATU for Catalysts and Solvents
        self.CatSolEncoder = nn.Sequential(
            self.GATU_CatSol,
            nn.Linear(self.GE_OutSize, int(emb_size/2))
        )

        # 1.2 FingerPrintEncoder Parser
        FP_params = ["FPSize", "Channels"]
        for key in FPEncoder.keys():
            if key not in FP_params:
                raise Exception("%s is not the param in FingerPrint Encoder")

            if key == "FPSize":
                FPSize = FPEncoder[key] # len_drfp

            if key == "Channels":
                FP_Channels = FPEncoder[key]  # [256, 256]
                assert len(FP_Channels) == 2

        # FPEncoder
        self.FPEncoder = nn.Sequential(
            nn.Linear(FPSize, FP_Channels[0]),
            nn.ReLU(),
            nn.Linear(FP_Channels[0], FP_Channels[1]),
            nn.ReLU(),
            nn.Linear(FP_Channels[1], emb_size)
        )

        # 3. Modality Alignment Parser
        MA_params = ["Heads"]
        for key in ModalityAlignment.keys():
            if key not in MA_params:
                raise Exception("%s is not the param in ModalityAlignment")

            if key == "Heads":
                MA_Heads = ModalityAlignment[key]  # 4

        self.MA = nn.TransformerEncoderLayer(d_model=emb_size * 2, nhead=MA_Heads)

        # 4. Decoder Parser
        Decoder_params = ["Channels"]
        for key in Decoder.keys():
            if key not in Decoder_params:
                raise Exception("%s is not the param in Decoder")

            if key == "Channels":
                Decoder_Channels = Decoder[key]  # [1000, 500, 100]
                assert len(Decoder_Channels) == 3

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(emb_size * 2, Decoder_Channels[0]),
            nn.ReLU(),
            nn.BatchNorm1d(Decoder_Channels[0]),

            nn.Linear(Decoder_Channels[0], Decoder_Channels[1]),
            nn.ReLU(),

            nn.Linear(Decoder_Channels[1], Decoder_Channels[2]),
            nn.ReLU(),

            nn.Linear(Decoder_Channels[2], 1)
        )


    def forward(self, x, y=None, epoch=None, device=None):
        # import data
        ReaPro_data, CatSol_data, FP = x
        # Graph Modality
        graph_emb = torch.cat([self.ReaProEncoder(ReaPro_data),
                                   self.CatSolEncoder(CatSol_data)], dim=1)
        # FP Modality
        FP_embed = self.FPEncoder(FP)
        x = torch.cat([graph_emb, FP_embed], dim=1)


        # Modality Alignment
        x = self.MA(x)

        # FDS operation
        if self.FDS is not None:
            if epoch >= 1:
                self.FDS.cpu()
                x = self.FDS.smooth(x.cpu(), y.cpu(), epoch).to(device)
                y = y.to(device)

        # Decoder
        x = self.decoder(x)

        return x

    def GetLatentVector(self, x):
        # import data
        ReaPro_data, CatSol_data, RxnSmi = x

        return self.ReaProEncoder(ReaPro_data), self.CatSolEncoder(CatSol_data), self.FPEncoder(RxnSmi)

# Model Evaluation Function
def RMSE(pred, true):
    return np.sqrt(mean_squared_error(true, pred))

def R2(pred, true):
   return r2_score(true, pred)

def MAE(pred, true):
    return mean_absolute_error(true, pred)