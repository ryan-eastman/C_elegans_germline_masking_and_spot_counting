"""The net3d network: a compact residual 3D U-Net that turns a (foci, DNA) volume into a per-voxel
heatmap of RAD-51 foci (one logit per voxel; the sigmoid is the probability that a focus centre is
there).

Design, for the 0.108 x 0.108 x 0.2 um grid of the lab's stacks (z voxels 1.85 times the xy ones):

* the first pooling is xy only (1 x 2 x 2), which leaves the next level at 0.2 x 0.217 x 0.217 um,
  roughly isotropic, and every later pooling is 2 x 2 x 2;
* every level is one residual block (two 3 x 3 x 3 convolutions with batch normalisation and SiLU, a
  1 x 1 x 1 projection on the skip when the width changes);
* pooling is a max over non-overlapping windows written as a reshape and ``amax``, and upsampling is a
  transposed convolution whose kernel equals its stride. Both have deterministic CUDA backward passes,
  which ``torch.nn.MaxPool3d`` and trilinear interpolation do not, so training can run under
  ``torch.use_deterministic_algorithms(True)``;
* batch normalisation, not group or instance normalisation: in eval mode it is a fixed per-channel
  affine map, so the output at a voxel never depends on how the volume was cut into tiles (the tiled
  inference of infer.py relies on that).

The output of a voxel depends only on the input inside its receptive field; ``receptive_radius``
gives a safe upper bound per axis, and the tiled inference uses it as the tile overlap so tiled and
untiled predictions agree to float rounding.
"""
from __future__ import annotations

import math

import torch
from torch import nn

DEFAULT_ARCH = {"in_channels": 2, "base_channels": 16, "levels": 4}


def pools_for(levels: int) -> list[tuple[int, int, int]]:
    """The pooling windows (z, y, x) between the ``levels`` levels: xy only first, then 2 x 2 x 2."""
    if int(levels) < 2:
        raise ValueError(f"the U-Net needs at least 2 levels, got {levels}")
    return [(1, 2, 2)] + [(2, 2, 2)] * (int(levels) - 2)


def size_divisor(levels: int) -> tuple[int, int, int]:
    """Every input extent (z, y, x) must be a multiple of this, and tile origins too (so the pooling
    windows of a tile coincide with those of the whole volume)."""
    d = [1, 1, 1]
    for p in pools_for(levels):
        d = [a * b for a, b in zip(d, p)]
    return tuple(d)


def receptive_radius(levels: int) -> tuple[int, int, int]:
    """An upper bound, per axis in voxels, on how far an input voxel can be from an output voxel it
    influences. Each 3 x 3 x 3 convolution adds one step of the current level on each side, each
    pooling and each transposed convolution with window k at most (k - 1) steps of the finer level
    (the window alignment can put the whole extra on one side). Exact tiling needs an overlap of at
    least this much."""
    pools = pools_for(levels)
    rad = [0, 0, 0]
    jump = [1, 1, 1]
    # encoder: a residual block (two convolutions) at every level, a pooling before every level but the first
    for lev in range(int(levels)):
        if lev > 0:
            p = pools[lev - 1]
            for a in range(3):
                rad[a] += (p[a] - 1) * jump[a]
                jump[a] *= p[a]
        for a in range(3):
            rad[a] += 2 * jump[a]
    # decoder: a transposed convolution then a residual block at every level but the deepest
    for lev in range(int(levels) - 1, 0, -1):
        p = pools[lev - 1]
        for a in range(3):
            jump[a] //= p[a]
            rad[a] += (p[a] - 1) * jump[a]
            rad[a] += 2 * jump[a]
    return tuple(int(r) for r in rad)


def tile_margin(levels: int) -> tuple[int, int, int]:
    """The tile overlap of the inference: the receptive radius rounded up to the size divisor."""
    div = size_divisor(levels)
    return tuple(int(math.ceil(r / d) * d) for r, d in zip(receptive_radius(levels), div))


def max_pool(x: torch.Tensor, window: tuple[int, int, int]) -> torch.Tensor:
    """Max over non-overlapping (kz, ky, kx) windows, as a reshape and an amax (deterministic forward
    and backward on every device; the extents must be multiples of the window)."""
    kz, ky, kx = window
    if (kz, ky, kx) == (1, 1, 1):
        return x
    n, c, d, h, w = x.shape
    if d % kz or h % ky or w % kx:
        raise ValueError(f"tensor extents {(d, h, w)} are not multiples of the pooling window {window}")
    return x.reshape(n, c, d // kz, kz, h // ky, ky, w // kx, kx).amax(dim=(3, 5, 7))


class ResBlock(nn.Module):
    """Two 3 x 3 x 3 convolutions with batch normalisation and SiLU, plus the (projected) input."""

    def __init__(self, cin: int, cout: int):
        super().__init__()
        self.conv1 = nn.Conv3d(cin, cout, 3, padding=1, bias=False)
        self.bn1 = nn.BatchNorm3d(cout)
        self.conv2 = nn.Conv3d(cout, cout, 3, padding=1, bias=False)
        self.bn2 = nn.BatchNorm3d(cout)
        self.skip = (nn.Identity() if cin == cout else
                     nn.Sequential(nn.Conv3d(cin, cout, 1, bias=False), nn.BatchNorm3d(cout)))
        self.act = nn.SiLU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.act(self.bn1(self.conv1(x)))
        h = self.bn2(self.conv2(h))
        return self.act(h + self.skip(x))


class Net3D(nn.Module):
    """The residual 3D U-Net. Input (N, in_channels, Z, Y, X) with Z, Y, X multiples of
    ``size_divisor(levels)``; output (N, 1, Z, Y, X) logits."""

    def __init__(self, in_channels: int = 2, base_channels: int = 16, levels: int = 4):
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
        self.head = nn.Conv3d(chans[0], 1, 1)
        # a focus centre is rare: start every voxel near p = 0.018 so the first steps do not have to
        # unlearn a 0.5 prior on half a million background voxels per patch
        nn.init.constant_(self.head.bias, -4.0)

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


def build_model(arch: dict | None = None) -> Net3D:
    a = {**DEFAULT_ARCH, **(arch or {})}
    return Net3D(in_channels=int(a["in_channels"]), base_channels=int(a["base_channels"]),
                 levels=int(a["levels"]))


def n_parameters(model: nn.Module) -> int:
    return int(sum(p.numel() for p in model.parameters()))
