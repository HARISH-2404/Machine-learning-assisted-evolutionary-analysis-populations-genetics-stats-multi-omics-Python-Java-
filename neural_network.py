"""Ten PyTorch deep-learning architectures and a basic trainer."""
def build_model(name,input_dim,output_dim=1,hidden_dim=128,n_layers=2,dropout=.2):
    try: import torch
    except ImportError as e: raise ImportError("pip install torch") from e
    nn=torch.nn
    class MLP(nn.Module):
        def __init__(self):
            super().__init__(); layers=[]; d=input_dim
            for _ in range(n_layers): layers += [nn.Linear(d,hidden_dim),nn.ReLU(),nn.Dropout(dropout)]; d=hidden_dim
            layers += [nn.Linear(d,output_dim)]; self.net=nn.Sequential(*layers)
        def forward(self,x): return self.net(x)
    class Residual(nn.Module):
        def __init__(self):
            super().__init__(); self.i=nn.Linear(input_dim,hidden_dim); self.blocks=nn.Sequential(*[nn.Sequential(nn.Linear(hidden_dim,hidden_dim),nn.ReLU(),nn.Linear(hidden_dim,hidden_dim)) for _ in range(n_layers)]); self.o=nn.Linear(hidden_dim,output_dim)
        def forward(self,x):
            h=torch.relu(self.i(x))
            for b in self.blocks: h=torch.relu(h+b(h))
            return self.o(h)
    class Conv(nn.Module):
        def __init__(self):
            super().__init__(); self.c=nn.Sequential(nn.Conv1d(1,32,5,padding=2),nn.ReLU(),nn.Conv1d(32,64,3,padding=1),nn.ReLU(),nn.AdaptiveAvgPool1d(1)); self.o=nn.Linear(64,output_dim)
        def forward(self,x): return self.o(self.c(x.unsqueeze(1)).squeeze(-1))
    class RNN(nn.Module):
        def __init__(self,kind):
            super().__init__(); self.r=(nn.LSTM if kind=="lstm" else nn.GRU)(1,hidden_dim,batch_first=True); self.o=nn.Linear(hidden_dim,output_dim)
        def forward(self,x): return self.o(self.r(x.unsqueeze(-1))[0][:,-1,:])
    class Transformer(nn.Module):
        def __init__(self):
            super().__init__(); self.p=nn.Linear(1,hidden_dim); self.t=nn.TransformerEncoder(nn.TransformerEncoderLayer(hidden_dim,4,batch_first=True),n_layers); self.o=nn.Linear(hidden_dim,output_dim)
        def forward(self,x): return self.o(self.t(self.p(x.unsqueeze(-1))).mean(1))
    class AE(nn.Module):
        def __init__(self,variational=False):
            super().__init__(); self.v=variational; self.e=nn.Sequential(nn.Linear(input_dim,hidden_dim),nn.ReLU(),nn.Linear(hidden_dim,hidden_dim//2),nn.ReLU()); self.mu=nn.Linear(hidden_dim//2,hidden_dim//4); self.lv=nn.Linear(hidden_dim//2,hidden_dim//4); self.d=nn.Sequential(nn.Linear(hidden_dim//4,hidden_dim//2),nn.ReLU(),nn.Linear(hidden_dim//2,input_dim))
        def forward(self,x):
            h=self.e(x); z=self.mu(h)
            if self.v: z=z+torch.randn_like(z)*torch.exp(.5*self.lv(h))
            return self.d(z)
    class Multi(nn.Module):
        def __init__(self):
            super().__init__(); self.s=nn.Sequential(nn.Linear(input_dim,hidden_dim),nn.ReLU(),nn.Dropout(dropout)); self.h=nn.ModuleList([nn.Linear(hidden_dim,1) for _ in range(output_dim)])
        def forward(self,x): return torch.cat([h(self.s(x)) for h in self.h],1)
    table={"mlp":lambda:MLP(),"residual_mlp":lambda:Residual(),"cnn1d":lambda:Conv(),"lstm":lambda:RNN("lstm"),"gru":lambda:RNN("gru"),"transformer":lambda:Transformer(),"autoencoder":lambda:AE(),"variational_autoencoder":lambda:AE(True),"denoising_autoencoder":lambda:AE(),"multi_task_mlp":lambda:Multi()}
    if name not in table: raise ValueError(f"Choose from {list(table)}")
    return table[name]()
def train_torch_model(model,X,y,epochs=100,batch_size=32,learning_rate=.001,task="regression",seed=42):
    try: import torch
    except ImportError as e: raise ImportError("pip install torch") from e
    import numpy as np
    torch.manual_seed(seed); X=torch.tensor(np.asarray(X),dtype=torch.float32); y=torch.tensor(np.asarray(y),dtype=torch.float32)
    if y.ndim==1:y=y[:,None]
    loader=torch.utils.data.DataLoader(torch.utils.data.TensorDataset(X,y),batch_size=batch_size,shuffle=True)
    opt=torch.optim.Adam(model.parameters(),lr=learning_rate); lossfn=torch.nn.MSELoss() if task=="regression" else torch.nn.BCEWithLogitsLoss(); hist=[]
    for _ in range(epochs):
        model.train(); total=0
        for xb,yb in loader:
            opt.zero_grad(); loss=lossfn(model(xb),yb); loss.backward(); opt.step(); total+=loss.item()*len(xb)
        hist.append(total/len(X))
    return {"model":model,"training_loss":hist}
