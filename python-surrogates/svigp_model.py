import numpy as np
import torch
import gpytorch
from torch.utils.data import TensorDataset, DataLoader

torch.set_default_dtype(torch.float32)

class SVGPModel(gpytorch.models.ApproximateGP):
    def __init__(self, inducing_points):
        variational_distribution = gpytorch.variational.CholeskyVariationalDistribution(
            inducing_points.size(0)
        )
        variational_strategy = gpytorch.variational.VariationalStrategy(
            self,
            inducing_points,
            variational_distribution,
            learn_inducing_locations=True
        )
        super().__init__(variational_strategy)

        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(
            gpytorch.kernels.RBFKernel(ard_num_dims=inducing_points.shape[1])
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def train_svgp(X_np,
               y_np,
               num_inducing=100,
               epochs=300,
               batch_size=256,
               lr=0.01,
               seed=123,
               device='cpu'):

    np.random.seed(seed)
    torch.manual_seed(seed)

    X_np = np.asarray(X_np, dtype=np.float32)
    y_np = np.asarray(y_np, dtype=np.float32).reshape(-1)

    n, d = X_np.shape
    num_inducing = int(min(num_inducing, n))
    batch_size = int(min(batch_size, n))

    # Standardize inputs and outputs
    xm = X_np.mean(axis=0)
    xs = X_np.std(axis=0)
    xs[xs == 0] = 1.0

    ym = float(y_np.mean())
    ys = float(y_np.std())
    if ys == 0:
        ys = 1.0

    Xs = (X_np - xm) / xs
    ys_std = (y_np - ym) / ys

    X = torch.tensor(Xs, dtype=torch.float32, device=device)
    y = torch.tensor(ys_std, dtype=torch.float32, device=device)

    # Simple inducing point initialization: random subset
    ind_idx = np.random.choice(n, size=num_inducing, replace=False)
    inducing_points = X[ind_idx, :].clone()

    model = SVGPModel(inducing_points).to(device)
    likelihood = gpytorch.likelihoods.GaussianLikelihood().to(device)

    model.train()
    likelihood.train()

    dataset = TensorDataset(X, y)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

    optimizer = torch.optim.Adam(
        [
            {'params': model.parameters()},
            {'params': likelihood.parameters()},
        ],
        lr=lr
    )

    mll = gpytorch.mlls.VariationalELBO(
        likelihood,
        model,
        num_data=n
    )

    for epoch in range(int(epochs)):
        for xb, yb in loader:
            optimizer.zero_grad()
            output = model(xb)
            loss = -mll(output, yb)
            loss.backward()
            optimizer.step()

    return {
        'model': model,
        'likelihood': likelihood,
        'xm': xm,
        'xs': xs,
        'ym': ym,
        'ys': ys
    }


def sample_svgp(obj,
                X_np,
                B=1000,
                device='cpu'):

    model = obj['model']
    likelihood = obj['likelihood']
    xm = obj['xm']
    xs = obj['xs']
    ym = obj['ym']
    ys = obj['ys']

    X_np = np.asarray(X_np, dtype=np.float32)
    Xs = (X_np - xm) / xs

    Xt = torch.tensor(Xs, dtype=torch.float32, device=device)

    model.eval()
    likelihood.eval()

    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        pred_dist = likelihood(model(Xt))
        samp = pred_dist.rsample(sample_shape=torch.Size([int(B)]))

    # Un-standardize
    samp = samp.cpu().numpy() * ys + ym
    return samp
