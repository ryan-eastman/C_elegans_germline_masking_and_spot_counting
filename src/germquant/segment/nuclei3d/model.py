"""The nuclei3d network: a compact residual 3D U-Net that turns a DAPI volume into three per-voxel maps.

Output channels (one 1 x 1 x 1 convolution on the last decoder level):

0. foreground logit: the sigmoid is the probability that the voxel lies inside a nucleus;
1. cell distance in micrometres: the distance from the voxel to the boundary of its own nucleus (0 in
   the background), the "cell distance" of Scherr et al. 2020 kept in physical units;
2. neighbour closeness, 0 to 1: 1 where another nucleus touches, falling to 0 at the neighbour radius
   (the "neighbour distance" of Scherr et al. 2020, inverted and steepened); 0 in the background.

The encoder and decoder are those of the net3d foci detector (germquant.spots.net3d.model, reused here,
not copied): four levels of one residual block each (two 3 x 3 x 3 convolutions with batch
normalisation and SiLU), an xy-only first pooling (1 x 2 x 2) for the 1.85 anisotropy of the 0.108 x
0.108 x 0.2 um grid, 2 x 2 x 2 pooling after it, pooling written as a reshape and a max and upsampling as
a transposed convolution whose kernel equals its stride (both deterministic on CUDA), and batch
normalisation (a fixed affine map in eval mode, so a voxel's output never depends on how the volume was
tiled). ``receptive_radius`` and ``tile_margin`` are the net3d functions: a 1 x 1 x 1 head adds nothing.
"""
from __future__ import annotations

import torch
from torch import nn

from germquant.spots.net3d.model import (  # noqa: F401  (re-exported for the inference and the tests)
    ResBlock,
    max_pool,
    pools_for,
    receptive_radius,
    size_divisor,
    tile_margin,
)

HEADS = ("foreground", "cell_distance_um", "neighbour")
DEFAULT_ARCH = {"in_channels": 1, "base_channels": 16, "levels": 4}
# initial head biases: about 12 percent of a training patch is foreground, so start the foreground
# logit near that prior; the two regression heads start at 0 (the background value)
HEAD_BIAS = (-2.0, 0.0, 0.0)


class Nuclei3DNet(nn.Module):
    """The residual 3D U-Net. Input (N, in_channels, Z, Y, X) with Z, Y, X multiples of
    ``size_divisor(levels)``; output (N, 3, Z, Y, X): foreground logit, cell distance (um), neighbour
    closeness (raw regression outputs; the inference clamps the two distances to their ranges)."""

    def __init__(self, in_channels: int = 1, base_channels: int = 16, levels: int = 4):
        super().__init__()
        self.in_channels = int(in_channels)
        self.base_channels = int(base_channels)
        self.levels = int(levels)
        self.pools = pools_for(self.levels)
        chans = [self.base_channels * 2 ** i for i in range(self.levels)]
        self.enc = nn.ModuleList([ResBlock(self.in_channels, chans[0])] +
                                 [ResBlock(chans[i - 1], chans[i]) for i in range(1, self.levels)])
        self.up = nn.ModuleList()
        self.dec = nn.ModuleList()
        for i in range(self.levels - 1, 0, -1):
            p = self.pools[i - 1]
            self.up.append(nn.ConvTranspose3d(chans[i], chans[i - 1], kernel_size=p, stride=p))
            self.dec.append(ResBlock(2 * chans[i - 1], chans[i - 1]))
        self.head = nn.Conv3d(chans[0], len(HEADS), 1)
        with torch.no_grad():
            self.head.bias.copy_(torch.tensor(HEAD_BIAS, dtype=self.head.bias.dtype))

    def arch(self) -> dict:
        return {"in_channels": self.in_channels, "base_channels": self.base_channels, "levels": self.levels}

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        skips = []
        h = x
        for i, block in enumerate(self.enc):
            if i > 0:
                h = max_pool(h, self.pools[i - 1])
            h = block(h)
            skips.append(h)
        for j, (up, dec) in enumerate(zip(self.up, self.dec)):
            h = up(h)
            h = dec(torch.cat([h, skips[-2 - j]], dim=1))
        return self.head(h)


def build_model(arch: dict | None = None) -> Nuclei3DNet:
    a = {**DEFAULT_ARCH, **(arch or {})}
    return Nuclei3DNet(in_channels=int(a["in_channels"]), base_channels=int(a["base_channels"]),
                       levels=int(a["levels"]))


def n_parameters(model: nn.Module) -> int:
    return int(sum(p.numel() for p in model.parameters()))
