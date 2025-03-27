# -*- coding: utf-8 -*-
"""
Author: Philipp Seidl
        ELLIS Unit Linz, LIT AI Lab, Institute for Machine Learning
        Johannes Kepler University Linz
Contact: seidl@ml.jku.at

Model related functionality
"""
from utils import top_k_accuracy
from plotutils import plot_loss, plot_topk, plot_nte
from molutils import convert_smiles_to_fp
import os
import numpy as np
import torch
import torch.nn as nn
from collections import defaultdict
from scipy import sparse
import logging
from tqdm import tqdm

log = logging.getLogger(__name__)

class ChemRXNDataset(torch.utils.data.Dataset):
    "Torch Dataset for ChemRXN containing Xs: the input as np array, target: the target molecules (or nothing), and ys: the label"
    def __init__(self, Xs, target, ys, is_smiles=False, fp_size=2048, fingerprint_type='morgan'):
        self.is_smiles=is_smiles
        if is_smiles:
            self.Xs = Xs
            self.target = target
            self.fp_size = fp_size
            self.fingerprint_type = fingerprint_type
        else:
            self.Xs = Xs.astype(np.float32)
            self.target = target.astype(np.float32)
        self.ys = ys
        self.ys_is_sparse = isinstance(self.ys, sparse.csr.csr_matrix)

    def __getitem__(self, k):
        mol_fp = self.Xs[k]
        if self.is_smiles:
            mol_fp = convert_smiles_to_fp(mol_fp, fp_size=self.fp_size, which=self.fingerprint_type).astype(np.float32)
        
        target = None if self.target is None else self.target[k]
        if self.is_smiles and self.target:
            target = convert_smiles_to_fp(target, fp_size=self.fp_size, which=self.fingerprint_type).astype(np.float32)
        
        label = self.ys[k]
        if isinstance(self.ys, sparse.csr.csr_matrix):
            label = label.toarray()[0]
            
        return (mol_fp, target, label)

    def __len__(self):
        return len(self.Xs)

class ModelConfig(object):
    def __init__(self, **kwargs):
        self.fingerprint_type = kwargs.pop("fingerprint_type", 'morgan')
        self.template_fp_type = kwargs.pop("template_fp_type", 'rdk')
        self.num_templates = kwargs.pop("num_templates", 401)
        self.fp_size = kwargs.pop("fp_size", 2048)
        self.fp_radius = kwargs.pop("fp_radius", 4)

        self.device = kwargs.pop("device", 'cuda' if torch.cuda.is_available() else 'cpu')
        self.batch_size = kwargs.pop("batch_size", 32)
        self.pooling_operation_state_embedding = kwargs.pop('pooling_operation_state_embedding', 'mean')
        self.pooling_operation_head = kwargs.pop('pooling_operation_head', 'max')

        self.dropout = kwargs.pop('dropout', 0.0)

        self.lr = kwargs.pop('lr', 1e-4)
        self.optimizer = kwargs.pop("optimizer", "Adam")

        self.activation_function = kwargs.pop('activation_function', 'ReLU')
        self.verbose = kwargs.pop("verbose", False)  # debugging or printing additional warnings / information set to True

        self.hopf_input_size = kwargs.pop('hopf_input_size', 2048)
        self.hopf_output_size = kwargs.pop("hopf_output_size", 768)
        self.hopf_num_heads = kwargs.pop("hopf_num_heads", 1)
        self.hopf_asso_dim = kwargs.pop("hopf_asso_dim", 768)
        self.hopf_association_activation = kwargs.pop("hopf_association_activation", None)
        self.hopf_beta = kwargs.pop("hopf_beta", 0.125)  #  1/(self.hopf_asso_dim**(1/2) sqrt(d_k)
        self.norm_input = kwargs.pop("norm_input", False)
        self.norm_asso = kwargs.pop("norm_asso", False)

        # additional experimental hyperparams
        if 'hopf_n_layers' in kwargs.keys():
            self.hopf_n_layers = kwargs.pop('hopf_n_layers', 0)
        if 'mol_encoder_layers' in kwargs.keys():
            self.mol_encoder_layers = kwargs.pop('mol_encoder_layers', 1)
        if 'temp_encoder_layers' in kwargs.keys():
            self.temp_encoder_layers = kwargs.pop('temp_encoder_layers', 1)
        if 'encoder_af' in kwargs.keys():
            self.encoder_af = kwargs.pop('encoder_af', 'ReLU')

        # additional kwargs
        for key, value in kwargs.items():
            try:
                setattr(self, key, value)
            except AttributeError as err:
                log.error(f"Can't set {key} with value {value} for {self}")
                raise err


class Encoder(nn.Module):
    """Simple FFNN with Layer Normalization and Residual Connections"""
    def __init__(self, input_size: int = 2048, output_size: int = 1024,
                    num_layers: int = 1, dropout: float = 0.3, af_name: str = 'None', 
                    norm_in: bool = False, norm_out: bool = False):
        super().__init__()
        self.ws = []
        self.setup_af(af_name)
        self.norm_in = (lambda k: k) if not norm_in else nn.LayerNorm(input_size, elementwise_affine=False)
        self.norm_out = (lambda k: k) if not norm_out else nn.LayerNorm(output_size, elementwise_affine=False)
        self.setup_ff(input_size, output_size, num_layers)
        self.dropout = nn.Dropout(p=dropout)

    def forward(self, x: torch.Tensor):
        residual = x
        x = self.norm_in(x)
        for i, w in enumerate(self.ws):
            if i == (len(self.ws) - 1):
                x = self.dropout(w(x))  # all except last have activation function
            else:
                x = self.dropout(self.af(w(x)))
        x = self.norm_out(x)
        x += residual  # Residual connection
        return x

    def setup_ff(self, input_size: int, output_size: int, num_layers=1):
        """Setup feed-forward NN with n-layers"""
        for n in range(num_layers):
            w = nn.Linear(input_size if n == 0 else output_size, output_size)
            torch.nn.init.kaiming_normal_(w.weight, mode='fan_in', nonlinearity='linear')  # equivalent to LeCun init
            setattr(self, f'W_{n}', w)
            self.ws.append(getattr(self, f'W_{n}'))

    def setup_af(self, af_name: str):
        """Set activation function"""
        if af_name is None or af_name == 'None':
            self.af = lambda k: k
        else:
            try:
                self.af = getattr(nn, af_name)()
            except AttributeError as err:
                log.error(f"Can't find activation-function {af_name} in torch.nn")
                raise err


class MoleculeEncoder(Encoder):
    """
    Class for Molecule encoder: can be any class mapping Smiles to a Vector (preferable differentiable ;)
    """
    def __init__(self, config):
        self.config = config


class FPMolEncoder(Encoder):
    """
    Fingerprint Based Molecular encoder
    """
    def __init__(self, config):
        super().__init__(input_size=config.hopf_input_size * config.hopf_num_heads,
                         output_size=config.hopf_asso_dim * config.hopf_num_heads,
                         num_layers=config.mol_encoder_layers,
                         dropout=config.dropout,
                         af_name=config.encoder_af,
                         norm_in=config.norm_input,
                         norm_out=config.norm_asso)
        self.config = config

    def forward_smiles(self, list_of_smiles: list):
        fp_tensor = self.convert_smiles_to_tensor(list_of_smiles)
        return self.forward(fp_tensor)

    def convert_smiles_to_tensor(self, list_of_smiles):
        fps = convert_smiles_to_fp(list_of_smiles, fp_size=self.config.fp_size, 
                                   which=self.config.fingerprint_type, radius=self.config.fp_radius)
        fps_tensor = torch.from_numpy(fps.astype(np.float32)).to(dtype=torch.float).to(self.config.device)
        return fps_tensor


class TemplateEncoder(Encoder):
    """
    Class for Template encoder: can be any class mapping a Smarts-Reaction to a Vector (preferable differentiable ;)
    """
    def __init__(self, config):
        super().__init__(input_size=config.hopf_input_size * config.hopf_num_heads,
                         output_size=config.hopf_asso_dim * config.hopf_num_heads,
                         num_layers=config.temp_encoder_layers,
                         dropout=config.dropout,
                         af_name=config.encoder_af,
                         norm_in=config.norm_input,
                         norm_out=config.norm_asso)
        self.config = config
        if config.temp_encoder_layers == 0:
            raise ValueError(f"template_encoder must have at least 1 layer")


class MHN(nn.Module):
    """
    Class for MHN Module: setup with Hopfield-pooling (Hopfield nets can be seen as an associative memory)

    Params:
    ---------------
    use_template_encoder: bool
        Flag to indicate if template encoder should be used or not
    layer2weight: bool
        Flag to indicate if weights should be layer-wise or globally shared
    """
    def __init__(self, config, use_template_encoder=True, layer2weight=True):
        super().__init__()
        self.config = config
        self.use_template_encoder = use_template_encoder
        self.layer2weight = layer2weight

        self.molecule_encoder = FPMolEncoder(config=config)
        self.template_encoder = TemplateEncoder(config=config)

        self.layernorm1 = nn.LayerNorm(config.hopf_input_size * config.hopf_num_heads)
        self.layernorm2 = nn.LayerNorm(config.hopf_asso_dim * config.hopf_num_heads)

    def set_templates(self, template_list, which='rdk', fp_size=None, radius=2, learnable=False, njobs=1, only_templates_in_batch=False):
        self.which = which
        self.fp_size = fp_size if fp_size else self.config.fp_size
        self.radius = radius
        self.learnable = learnable
        self.only_templates_in_batch = only_templates_in_batch
        self.template_list = template_list
        self.config.num_templates = len(template_list)

        templates = convert_smiles_to_fp(template_list, fp_size=self.fp_size, radius=self.radius, njobs=njobs)
        templates_tensor = torch.from_numpy(templates.astype(np.float32)).to(dtype=torch.float).to(self.config.device)
        self.register_buffer("templates_tensor", templates_tensor, persistent=False)
        self.templates_tensor.requires_grad = learnable

    def set_templates_recursively(self):
        for i in range(self.config.hopf_n_layers):
            self.set_templates()

    def update_template_embedding(self, fp_size=2048, radius=4, which='rdk', learnable=False, njobs=1, only_templates_in_batch=False, template_list=None, verbose=True):
        if verbose:
            print('Updating Template embeddings...')

        if not template_list:
            template_list = self.template_list
        self.set_templates(template_list, which, fp_size, radius, learnable, njobs, only_templates_in_batch)

        return None

    def np_fp_to_tensor(self, np_fp):
        return torch.from_numpy(np_fp.astype(np.float64)).to(self.config.device).float()

    def masked_loss_fun(self, loss_fun, h_out, ys_batch):
        mask = torch.isfinite(ys_batch)
        if mask.sum() == 0:
            return torch.tensor(0.0, requires_grad=True)
        return loss_fun(h_out[mask], ys_batch[mask])

    def compute_losses(self, out, ys_batch, head_loss_weight=None):
        loss_fun = nn.CrossEntropyLoss()
        losses = defaultdict(list)
        for h_out, ys in zip(out, ys_batch):
            for i, h in enumerate(h_out):
                loss = self.masked_loss_fun(loss_fun, h, ys)
                if head_loss_weight:
                    loss *= head_loss_weight[i]
                losses[i].append(loss)

        return losses

    def forward_smiles(self, list_of_smiles, templates=None):
        m_tensor = self.molecule_encoder.forward_smiles(list_of_smiles)
        if templates is not None:
            templates = templates.to(self.config.device)
        return self.forward(m_tensor, templates)

    def encode_templates(self, list_of_smarts, batch_size=32, njobs=1):
        template_tensors = []
        for i in range(0, len(list_of_smarts), batch_size):
            batch_smarts = list_of_smarts[i:i + batch_size]
            batch_tensors = self.template_encoder.forward_smiles(batch_smarts)
            template_tensors.append(batch_tensors)

        return torch.cat(template_tensors, dim=0)

    def encode_smiles(self, list_of_smiles, batch_size=32, njobs=1):
        smiles_tensors = []
        for i in range(0, len(list_of_smiles), batch_size):
            batch_smiles = list_of_smiles[i:i + batch_size]
            batch_tensors = self.molecule_encoder.forward_smiles(batch_smiles)
            smiles_tensors.append(batch_tensors)

        return torch.cat(smiles_tensors, dim=0)

    def forward(self, m, templates=None):
        bs = m.shape[0]  # batch_size

        if (templates is None) and (self.templates is None):
            raise Exception('Either pass in templates, or init templates by running clf.set_templates')

        n_temp = len(templates) if templates is not None else len(self.templates)

        if self.training or (templates is not None) or (self.templates is None):
            templates = templates if templates is not None else self.templates
            X = self.template_encoder(templates)
        else:
            X = self.templates  # precomputed from last forward run

        Xi = self.mol_encoder(m)

        Xi = Xi.view(bs, self.config.hopf_num_heads, self.config.hopf_asso_dim)  # [bs, H, A]
        X = X.view(1, n_temp, self.config.hopf_asso_dim, self.config.hopf_num_heads)  # [1, T, A, H]

        XXi = torch.tensordot(Xi, X, dims=[(2, 1), (2, 0)])  # AxA -> [bs, T, H]

        # pooling over heads
        if self.config.hopf_num_heads <= 1:
            XXi = XXi[:, :, 0]  # torch.squeeze(QKt, dim=2)
        else:
            XXi = self.pooling_operation_head(XXi, dim=2)  # default is max pooling over H [bs, T]
            if (self.config.pooling_operation_head == 'max') or (self.config.pooling_operation_head == 'min'):
                XXi = XXi[0]  # max and min also return the indices =S

        out = self.beta * XXi  # [bs, T, H] # softmax over dim=1 #pooling_operation_head

        xinew = self.softmax(out) @ X.view(n_temp, self.config.hopf_asso_dim)  # [bs,T]@[T,emb] -> [bs,emb]

        if self.W_v:
            # call layers recursive
            hopfout = self.W_v(xinew)  # [bs,emb]@[emb,hopf_inp]  --> [bs, hopf_inp]
            # Residual connection
            hopfout = hopfout + m
            # Layer normalization
            hopfout = self.layernorm1(hopfout)
            # give it to the next layer
            out2 = self.layer.forward(hopfout)  # templates=self.W_v(self.K)
            # Residual connection
            out2 += m  # Updated here to add x instead of out
            # Layer normalization
            out2 = self.layernorm2(out2)
            # Combine with original output
            out = out * (1 - self.layer2weight) + out2 * self.layer2weight

        return out

    def train_from_np(self, Xs, targets, ys, is_smiles=False, epochs=2, lr=0.001, bs=32,
                      permute_batches=False, shuffle=True, optimizer=None,
                      use_dataloader=True, verbose=False,
                      wandb=None, scheduler=None, only_templates_in_batch=False):
        self.to(self.config.device)

        # Ensure templates_tensor is initialized before training
        if not hasattr(self, "templates_tensor") or self.templates_tensor is None:
            if self.template_list is None:
                raise ValueError("Template list not set. Use set_templates() before training.")
            self.update_template_embedding()

        dataset = ChemRXNDataset(Xs, targets, ys, is_smiles=is_smiles, fp_size=self.config.fp_size,
                                 fingerprint_type=self.config.fingerprint_type)
        if use_dataloader:
            dataloader = torch.utils.data.DataLoader(dataset, batch_size=bs, shuffle=shuffle)
        else:
            dataloader = dataset

        if optimizer is None:
            optimizer = torch.optim.Adam(self.parameters(), lr=lr)

        for epoch in range(epochs):
            self.train()
            epoch_loss = 0
            with tqdm(dataloader, unit="batch") as tepoch:
                for batch in tepoch:
                    tepoch.set_description(f"Epoch {epoch + 1}")

                    X_batch, target_batch, y_batch = batch
                    X_batch, y_batch = X_batch.to(self.config.device), y_batch.to(self.config.device)

                    optimizer.zero_grad()
                    outputs = self.forward(X_batch, templates=None)
                    loss_dict = self.compute_losses(outputs, y_batch)
                    loss = sum([l.mean() for l in loss_dict.values()])
                    loss.backward()
                    optimizer.step()
                    epoch_loss += loss.item()
                    tepoch.set_postfix(loss=epoch_loss)

                    if wandb:
                        wandb.log({"loss": epoch_loss})
            if scheduler:
                scheduler.step()

        return self

    def evaluate(self, Xs, targets, ys, split='test', is_smiles=False, bs=32, shuffle=False, wandb=None, only_loss=False):
        self.eval()

        dataset = ChemRXNDataset(Xs, targets, ys, is_smiles=is_smiles, fp_size=self.config.fp_size,
                                 fingerprint_type=self.config.fingerprint_type)
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=bs, shuffle=shuffle)

        total_loss = 0
        total_topk_accuracy = defaultdict(list)

        with torch.no_grad():
            for X_batch, target_batch, y_batch in dataloader:
                X_batch, y_batch = X_batch.to(self.config.device), y_batch.to(self.config.device)

                outputs = self.forward(X_batch, templates=None)
                loss_dict = self.compute_losses(outputs, y_batch)
                loss = sum([l.mean() for l in loss_dict.values()])
                total_loss += loss.item()

                if not only_loss:
                    for k in [1, 3, 5, 10, 20, 50, 100]:
                        topk_acc = top_k_accuracy(outputs, y_batch, k=k)
                        total_topk_accuracy[k].append(topk_acc)

        avg_loss = total_loss / len(dataloader)
        avg_topk_accuracy = {k: np.mean(v) for k, v in total_topk_accuracy.items()}

        if wandb:
            wandb.log({f"{split}_loss": avg_loss, **{f"{split}_top{k}_accuracy": v for k, v in avg_topk_accuracy.items()}})

        return avg_loss, avg_topk_accuracy

    def save_hist(self, results, output_dir=None, prefix=''):
        if output_dir is not None:
            os.makedirs(output_dir, exist_ok=True)
            result_file = os.path.join(output_dir, f"{prefix}results.pkl")
            print(f"Writing results to {result_file}")
            with open(result_file, 'wb') as f:
                torch.save(results, f)
        return

    def save_model(self, output_dir):
        if output_dir is not None:
            os.makedirs(output_dir, exist_ok=True)
            print(f"Saving model to {output_dir}")
            torch.save(self.state_dict(), os.path.join(output_dir, 'model.pt'))
            torch.save(self.config, os.path.join(output_dir, 'config.pt'))
        return


class SeglerBaseline(MHN):
    """FFNN - only the Molecule Encoder + an output projection"""
    def __init__(self, config=None):
        config.template_fp_type = 'none'
        config.temp_encoder_layers = 0
        super().__init__(config, use_template_encoder=False)
        self.W_out = torch.nn.Linear(config.hopf_asso_dim, config.num_templates)
        self.optimizer = getattr(torch.optim, self.config.optimizer)(self.parameters(), lr=self.lr)
        self.steps = 0
        self.hist = defaultdict(list)
        self.to(self.config.device)

    def forward(self, m, templates=None):
        """
        m: molecule in the form batch x fingerprint
        templates: won't be used in this case
        returns logits ranking the templates for each molecule
        """
        bs = m.shape[0] #batch_size
        Xi = self.mol_encoder(m)
        Xi = self.mol_encoder.af(Xi) # is not applied in encoder for last layer
        out = self.W_out(Xi) # [bs, T] # softmax over dim=1
        return out

class StaticQK(MHN):
    """ Static QK baseline - beware to have the same fingerprint for mol_encoder as for the template_encoder (fp2048 r4 rdk by default)"""
    def __init__(self, config=None):
        if config:
            self.config = config
        else:
            self.config = ModelConfig()
        super().__init__(config)

        self.fp_size = 2048
        self.fingerprint_type = 'rdk'
        self.beta = 1

    def update_template_embedding(self, which='rdk', fp_size=2048, radius=4, learnable=False):
        bs = self.config.batch_size
        split_template_list = [t.split('>>')[0].split('.') for t in self.template_list]
        self.templates = torch.from_numpy(convert_smiles_to_fp(split_template_list, 
                                                               is_smarts=True, fp_size=fp_size, 
                                                               radius=radius, which=which).max(1)).float().to(self.config.device)


    def forward(self, m, templates=None):
        """

        """
        #states_emb = self.fcfe(state_fp)
        bs = m.shape[0] #batch_size
        
        Xi = m #[bs, emb]
        X = self.templates #[T, emb])

        XXi = Xi@X.T # [bs, T]
        
        # normalize
        t_sum = templates.sum(1) #[T]
        t_sum = t_sum.view(1,-1).expand(bs, -1) #[bs, T]
        XXi = XXi / t_sum
        
        # not neccecaire because it is not trained
        out = self.beta*XXi # [bs, T] # softmax over dim=1
        return out

class Retrosim(StaticQK):
    """ Retrosim-like baseline only for template relevance prediction """
    def fit_with_train(self, X_fp_train, y_train):
        self.templates = torch.from_numpy(X_fp_train).float().to(self.config.device)
        # train_samples, num_templates
        self.sample2acttemplate = torch.nn.functional.one_hot(torch.from_numpy(y_train), self.config.num_templates).float()
        tmpnorm = self.sample2acttemplate.sum(0)
        tmpnorm[tmpnorm==0] = 1
        self.sample2acttemplate = (self.sample2acttemplate / tmpnorm).to(self.config.device) # results in an average after dot product

    def forward(self, m, templates=None):
        """
        """
        out = super().forward(m, templates=templates)
        # bs, train_samples

        # map out to actual templates
        out = out @ self.sample2acttemplate

        return out