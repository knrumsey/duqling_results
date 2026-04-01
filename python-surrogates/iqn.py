"""Implicit Quantile Network (IQN) — the GBC backbone.

Architecture: cosine quantile embedding (Dabney et al. 2018).
Loss: three-term (L1 mean anchor + monotonicity + quantile).
Optimizer: Adam + cosine annealing.

Code copied from: https://github.com/VadimSokolov/gbc-surrogate/
"""

import numpy as np
import torch
import torch.nn as nn


class IQN(nn.Module):
    """Implicit Quantile Network.

    Cosine embedding: phi(tau) = [cos(0*pi*tau), ..., cos((nh-1)*pi*tau)]
    Output: (batch, 2) where col0=mean_hat, col1=quantile_hat.
    """

    def __init__(self, xdim, hdim=256, nh=32):
        super().__init__()
        self.nh = nh
        self.fc_tau = nn.Sequential(nn.Linear(nh, hdim), nn.ReLU())
        self.fc_x = nn.Sequential(nn.Linear(xdim, hdim), nn.ReLU())
        self.fc1 = nn.Sequential(nn.Linear(hdim, hdim), nn.ReLU())
        self.fc2 = nn.Sequential(nn.Linear(hdim, 64), nn.Tanh())
        self.fc_out = nn.Linear(64, 2)
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x, tau):
        i = torch.arange(self.nh, dtype=torch.float32, device=x.device)
        h_tau = self.fc_tau(torch.cos(i * torch.pi * tau))
        h_x = self.fc_x(x)
        h = self.fc1(h_x * h_tau.unsqueeze(0))
        return self.fc_out(self.fc2(h))

    def loss_fn(self, x, y, w=(0.3, 0.3, 0.4)):
        """Three-term loss.

        w[0]: L1 on conditional mean (location anchor).
        w[1]: monotonicity regularization (quantile ordering).
        w[2]: standard quantile loss.
        """
        tau = torch.rand(1, device=x.device).item()
        tauind = tau < 0.5
        f = self(x, tau)
        e = y.view(-1, 1) - f
        loss = w[0] * torch.mean(torch.abs(e[:, 0]))
        mono = (tauind * torch.mean(torch.relu(-e[:, 1]))
                + (1 - tauind) * torch.mean(torch.relu(e[:, 1])))
        loss += w[1] * abs(tau - 0.5) * mono
        loss += w[2] * torch.mean(
            torch.maximum(tau * e[:, 1], (tau - 1) * e[:, 1]))
        return loss


def train_iqn(X_np, y_np, epochs=3000, hdim=256, nh=32, lr=1e-3, wd=1e-4,
              seed=42, loss_w=(0.3, 0.3, 0.4), device=None):
    """Train IQN with Adam + cosine annealing.

    Returns: (model, xm, xs, ym, ys) — model and normalization stats.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(seed)
    xm, xs = X_np.mean(0), X_np.std(0) + 1e-8
    ym, ys = float(y_np.mean()), float(y_np.std()) + 1e-8
    Xt = torch.tensor((X_np - xm) / xs, dtype=torch.float32, device=device)
    yt = torch.tensor((y_np - ym) / ys, dtype=torch.float32, device=device)

    model = IQN(X_np.shape[1], hdim=hdim, nh=nh).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(
        opt, T_max=epochs, eta_min=lr * 0.01)

    model.train()
    for _ in range(epochs):
        opt.zero_grad()
        loss = model.loss_fn(Xt, yt, w=loss_w)
        loss.backward()
        opt.step()
        sched.step()
    model.eval()
    return model, xm, xs, ym, ys


def sample_iqn(model, X_np, xm, xs, ym, ys, B=500, chunk=1000, device=None):
    """Sample B quantiles from trained IQN.

    Returns: (B, n) array of quantile predictions in original scale.
    """
    if device is None:
        device = next(model.parameters()).device
    Xt = torch.tensor((X_np - xm) / xs, dtype=torch.float32, device=device)
    taus = torch.linspace(0.005, 0.995, B, device=device)
    rows = []
    with torch.no_grad():
        for tau in taus:
            q_list = []
            for i in range(0, len(Xt), chunk):
                q_list.append(
                    model(Xt[i:i + chunk], tau.item())[:, 1].cpu().numpy())
            rows.append(np.concatenate(q_list) * ys + ym)
    return np.array(rows)
