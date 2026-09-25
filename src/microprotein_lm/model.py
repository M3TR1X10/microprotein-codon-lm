"""Small causal decoder, adapted from the ChaitanyaK77 SLM scaffold.

See references/SLM-LICENSE and docs/scaffold.md for attribution and changes.
"""
import math
from dataclasses import dataclass, asdict

import torch
from torch import nn
from torch.nn import functional as F


@dataclass
class ModelConfig:
    vocab_size: int
    max_tokens: int
    token_width: int
    d_model: int = 48
    n_heads: int = 4
    n_layers: int = 2
    dropout: float = 0.1

    def __post_init__(self):
        if self.vocab_size not in (4,64) or self.token_width not in (1,3):
            raise ValueError('Only fixed base and codon vocabularies are supported')
        if (self.vocab_size,self.token_width) not in ((4,1),(64,3)):
            raise ValueError('Vocabulary and token width disagree')
        if self.d_model % self.n_heads or self.d_model % 2 or self.d_model < 2:
            raise ValueError('d_model must be even and divisible by n_heads')
        if self.n_layers < 1 or self.max_tokens < 1 or not 0 <= self.dropout < 1:
            raise ValueError('Invalid model dimensions or dropout')


class Attention(nn.Module):
    def __init__(self, c):
        super().__init__()
        self.qkv = nn.Linear(c.d_model,3*c.d_model)
        self.proj = nn.Linear(c.d_model,c.d_model)
        self.heads = c.n_heads
        self.dropout = c.dropout
        self.residual_dropout = nn.Dropout(c.dropout)

    def forward(self,x):
        b,t,d = x.shape
        q,k,v = [z.view(b,t,self.heads,d//self.heads).transpose(1,2) for z in self.qkv(x).chunk(3,dim=-1)]
        y = F.scaled_dot_product_attention(q,k,v,is_causal=True,
            dropout_p=self.dropout if self.training else 0.0)
        return self.residual_dropout(self.proj(y.transpose(1,2).contiguous().view(b,t,d)))


class Block(nn.Module):
    def __init__(self,c):
        super().__init__()
        self.ln1 = nn.LayerNorm(c.d_model)
        self.ln2 = nn.LayerNorm(c.d_model)
        self.attn = Attention(c)
        self.mlp = nn.Sequential(nn.Linear(c.d_model,4*c.d_model),nn.GELU(),
            nn.Linear(4*c.d_model,c.d_model),nn.Dropout(c.dropout))

    def forward(self,x):
        x = x + self.attn(self.ln1(x))
        return x + self.mlp(self.ln2(x))


class Decoder(nn.Module):
    def __init__(self,c):
        super().__init__()
        self.config = c
        self.embedding = nn.Embedding(c.vocab_size,c.d_model)
        # Fixed positions measured in nucleotides; no extra learned positional
        # parameters in the longer base-token arm.
        position = torch.arange(c.max_tokens).unsqueeze(1)*c.token_width
        frequency = torch.exp(torch.arange(0,c.d_model,2)*(-math.log(10000)/c.d_model))
        pe = torch.zeros(c.max_tokens,c.d_model)
        pe[:,0::2] = torch.sin(position*frequency)
        pe[:,1::2] = torch.cos(position*frequency)
        self.register_buffer('position',pe)
        self.dropout = nn.Dropout(c.dropout)
        self.blocks = nn.ModuleList([Block(c) for _ in range(c.n_layers)])
        self.norm = nn.LayerNorm(c.d_model)
        self.head = nn.Linear(c.d_model,c.vocab_size,bias=False)
        self.head.weight = self.embedding.weight
        self.apply(self._init)

    def _init(self,module):
        if isinstance(module,(nn.Linear,nn.Embedding)):
            nn.init.normal_(module.weight,std=0.02)
        if isinstance(module,nn.Linear) and module.bias is not None:
            nn.init.zeros_(module.bias)

    def forward(self,x):
        if x.shape[1] > self.config.max_tokens:
            raise ValueError('Sequence exceeds configured biological context')
        x = self.dropout(self.embedding(x)+self.position[:x.shape[1]])
        for block in self.blocks:
            x = block(x)
        return self.head(self.norm(x))


def summed_loss(logits,targets):
    return F.cross_entropy(logits.reshape(-1,logits.shape[-1]),targets.reshape(-1),
        ignore_index=-100,reduction='sum')
