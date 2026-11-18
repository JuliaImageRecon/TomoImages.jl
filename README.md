# TomoImages.jl

TomoImages.jl describes where the voxels of tomographic images lie: in which coordinate system, at
which positions and in which orientation, and what their other dimensions are. Its types can be
used by reconstruction packages to return their images and by packages that register, segment or
combine images to read them, so that images from different scanners can be compared, resampled onto
each other and combined.

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaImageRecon.github.io/TomoImages.jl/dev)
[![Build Status](https://github.com/JuliaImageRecon/TomoImages.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaImageRecon/TomoImages.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaImageRecon/TomoImages.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaImageRecon/TomoImages.jl)

> [!IMPORTANT]
> TomoImages.jl is in an early stage of development. Its interface may change between versions.

- **Grid geometries** with an origin (the center of the first voxel), a spacing, a direction matrix
  for rotated and reflected grids, and the coordinate space they are written in.
- **Spaces** as types: a default space, the patient spaces LPS and RAS, and spaces that other
  packages define, for example for the coordinates of a scanner.
- **Images** that wrap any array (on the CPU or a GPU) with its geometry, a name for every
  dimension, coordinates or labels for the non-grid dimensions (time, contrast, coil) and metadata.
- **Transforms** between spaces, checked when they are composed or used.
- **Resampling** onto another grid, with or without a transform: nearest neighbor and multilinear
  interpolation as one KernelAbstractions kernel for CPUs and GPUs, B-splines with
  Interpolations.jl, and resamplers of other packages.

```julia
using TomoImages

g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=LPS())
img = TomoImage(rand(Float32, 64, 64, 32), g; dims=(x=GridAxis(1), y=GridAxis(2), z=GridAxis(3)))

index_to_world(g, (1, 1, 1))       # the origin: the center of the first voxel
slice = img[z=16]                  # a 64x64 image whose geometry is a slice in the 3D space

coarse = GridGeometry((32, 32, 16); spacing=(1.0, 1.0, 2.0), origin=(-15.5, -15.5, -15.0), space=LPS())
resample(img, coarse)              # img on the coarse grid, multilinear interpolation
```

The documentation (`docs/`, built with `julia --project=docs docs/make.jl`) has the background and
scope of the package, a getting-started page, tutorials for
reconstruction packages and for users of their images, how-to guides for spaces and resamplers, an
explanation of the conventions, the reference and a list of related packages.

## License

MIT, see [LICENSE](LICENSE).
