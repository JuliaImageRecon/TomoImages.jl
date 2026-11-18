# TomoImages.jl

TomoImages.jl describes where the voxels of tomographic images lie: in which coordinate system, at
which positions and in which orientation. Reconstruction packages return their images in its types,
so that images from different devices can be compared, resampled onto each other and combined.

- **Grid geometries** with an origin (the center of the first voxel), a spacing, a direction matrix
  for rotated and reflected grids, and the coordinate space they are written in.
- **Spaces** as types: device spaces, the patient spaces LPS and RAS, and spaces of other packages.
- **Images** that wrap any array (on the CPU or a GPU) with its geometry, named axes for the
  non-spatial dimensions (time, contrast, coil) and metadata.
- **Transforms** between spaces, checked when they are composed or used.
- **Resampling** onto another grid, with or without a transform: nearest neighbor and multilinear
  interpolation as one KernelAbstractions kernel for CPUs and GPUs, B-splines with
  Interpolations.jl, and resamplers of other packages.

```julia
using TomoImages

g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=LPS())
img = SpatialImage(rand(Float32, 64, 64, 32), g)

index_to_world(g, (1, 1, 1))       # the origin: the center of the first voxel
slice = img[:, :, 16]              # a 64x64 image whose geometry is a slice in the 3D space

coarse = GridGeometry((32, 32, 16); spacing=(1.0, 1.0, 2.0), origin=(-15.5, -15.5, -15.0), space=LPS())
resample(img, coarse)              # img on the coarse grid, multilinear interpolation
```

The documentation (`docs/`, built with `julia --project=docs docs/make.jl`) has a getting-started
page, a tutorial that overlays images from two devices, how-to guides for spaces and resamplers, an
explanation of the conventions, and the reference.

## License

MIT, see [LICENSE](LICENSE).
