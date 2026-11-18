# TomoImages.jl

TomoImages.jl describes where the voxels of tomographic images lie: in which coordinate system, at
which positions and in which orientation. Reconstruction packages return their images in its types,
so that images from different devices can be compared, resampled onto each other and combined.

- **Grid geometries** with an origin (the center of the first voxel), a spacing, a direction matrix
  for rotated and reflected grids, and the coordinate space they are written in. Any number of
  dimensions, and slices that stay in the space of their volume.
- **Spaces** as types: device spaces, the patient spaces LPS and RAS, and spaces of other packages.
  Images and transforms in different spaces are never combined without a transform.
- **Images** that wrap any array (on the CPU or a GPU) with its geometry, named axes for the
  non-spatial dimensions (time, contrast, coil), and metadata. Indexing, views and broadcasting keep
  the geometry.
- **Transforms** between spaces, checked when they are composed or used.
- **Resampling** onto another grid, with or without a transform: nearest neighbor and multilinear
  interpolation as one KernelAbstractions kernel for CPUs and GPUs, B-splines with Interpolations.jl,
  and resamplers of other packages.

## Example

```julia
using TomoImages

g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=LPS())
img = SpatialImage(rand(Float32, 64, 64, 32), g)

index_to_world(g, (1, 1, 1))       # the origin: the center of the first voxel
slice = img[:, :, 16]              # a 64x64 image whose geometry is a slice in the 3D space

coarse = GridGeometry((32, 32, 16); spacing=(1.0, 1.0, 2.0), origin=(-15.5, -15.5, -15.0), space=LPS())
resample(img, coarse)              # img on the coarse grid, multilinear interpolation
```

## Documentation

- [Getting Started](@ref) builds a grid and an image, selects parts of it and resamples it onto a
  rotated grid.
- The tutorial [Overlay Images from Two Devices](@ref) resamples an image of one device onto the
  grid of another through a transform.
- The how-to guides show how to [define a space](@ref "Define a Space") and how to
  [write a resampler](@ref "Write a Resampler").
- [Grids, Spaces and Transforms](@ref) explains the conventions: what the origin, spacing and
  direction mean, what a space is, and in which direction transforms are stored.
- The reference lists all types and functions: [Spaces and Grids](api/grids.md),
  [Transforms](api/transforms.md), [Images](api/images.md) and [Resampling](api/resampling.md).
