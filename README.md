# TomoImages.jl

This package describes where the voxels of tomographic images lie: in which coordinate system, at
which positions and in which orientation, and what their other dimensions are. Its types can be
used by reconstruction packages to return their images and by packages that register, segment or
combine images to read them, on the CPU and on GPUs.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaImageRecon.github.io/TomoImages.jl/dev)
[![Build Status](https://github.com/JuliaImageRecon/TomoImages.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaImageRecon/TomoImages.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaImageRecon/TomoImages.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaImageRecon/TomoImages.jl)

> [!IMPORTANT]
> TomoImages.jl is in an early stage of development. Its interface may change between versions.

## Features

- **Grid geometries**: an origin (the center of the first voxel), a spacing, a direction matrix for
  rotated and reflected grids, and the coordinate system (space) they are given in; any number of
  dimensions, and slices that keep their place in the space of their volume.
- **Spaces as types**: a default space, the patient spaces LPS and RAS, and spaces that other
  packages define, for example for the coordinates of a scanner. Images and transforms in different
  spaces are only combined through a transform that connects them.
- **Images with named dimensions**: a `TomoImage` wraps any array with its geometry, a name for
  every dimension, coordinates or labels for the dimensions that are not part of the grid (time,
  contrast, coil, echo), and metadata. Indexing, views, `permutedims`, `reverse` and broadcasting
  keep the description.
- **Transforms** between spaces, checked when they are composed or used.
- **Resampling** onto any other grid, through a transform where the spaces differ: nearest neighbor
  and multilinear interpolation as one
  [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) kernel for CPUs and
  GPUs, B-splines with [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl), and
  resamplers of other packages.
- **Custom components**: other packages add spaces, resamplers and image types through multiple
  dispatch, without changing TomoImages.
- Part of the [JuliaImageRecon](https://github.com/JuliaImageRecon) family, next to reconstruction
  packages such as [CTOperators.jl](https://github.com/JuliaImageRecon/CTOperators.jl).

## Example

```julia
using TomoImages

# 64 × 64 × 32 voxels of 0.5 × 0.5 × 1 mm in the patient space LPS; origin: the first voxel center
g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=LPS())
img = TomoImage(rand(Float32, 64, 64, 32), g; dims=(x=GridAxis(1), y=GridAxis(2), z=GridAxis(3)))

index_to_world(g, (1, 1, 1))       # the origin
slice = img[z=16]                  # a 64 × 64 image whose geometry is a slice in the 3D space

coarse = GridGeometry((32, 32, 16); spacing=(1.0, 1.0, 2.0), origin=(-15.5, -15.5, -15.0), space=LPS())
resample(img, coarse)              # img on the coarse grid, multilinear interpolation
```

On a GPU, wrap an array that lives on the GPU, e.g. a `CuArray`; resampling then runs there.

## Installation

TomoImages.jl is not registered yet:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaImageRecon/TomoImages.jl")
```

## Community Standards

This project is part of the Julia community and follows the
[Julia community standards](https://julialang.org/community/standards/).
