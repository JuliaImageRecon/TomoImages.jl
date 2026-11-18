# TomoImages.jl

TomoImages.jl describes reconstructed images together with their place in space: the positions
and orientation of their voxels, the coordinate system these are given in, and the meaning of every
dimension. Reconstruction packages return their images in its types, and packages that register,
segment or combine images read them, so that images from different scanners and modalities can be
compared and resampled onto each other, on the CPU and on GPUs.

```@example home
using TomoImages

g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=LPS())
img = SpatialImage(rand(Float32, 64, 64, 32), g; dims=(x=GridAxis(1), y=GridAxis(2), z=GridAxis(3)))
slice = img[z=16]                               # a 64 x 64 image that keeps its place in the volume
index_to_world(geometry(slice), (1, 1)), space(slice)
```

```@example home
coarse = GridGeometry((32, 32, 16); spacing=(1.0, 1.0, 2.0), origin=(-15.5, -15.5, -15.0), space=LPS())
summary(resample(img, coarse))                  # the image on a coarser grid, multilinear interpolation
```

!!! warning
    TomoImages.jl is in an early stage of development; its interface may still change.

## Background

A tomographic image is reconstructed from measurements: computed tomography (CT) from X-ray
projections, magnetic resonance imaging (MRI) from the signal of nuclear spins in k-space, magnetic
particle imaging (MPI) from the response of magnetic nanoparticles injected as a tracer, emission
tomography (PET, SPECT) from the decays of a radioactive tracer. Each reconstruction ends with an
array, and an array knows only indices. It does not record how far apart its voxels are, whether
the slab was tilted against the scanner, whether its positions refer to the coordinates of one
scanner, of another scanner or of the patient, or which of its dimensions are space and which are
time, contrast, coil or echo. Two images whose numbers look alike may cover different parts of the
body.

Every package that combines images, such as registration, segmentation, viewers or the overlay of
an MPI image on a CT image, needs this information and has to keep it consistent through every
slice, permutation and resampling. TomoImages keeps it with the image, much as packages for
geospatial data such as Rasters.jl keep the coordinate system and placement of maps (see
[Related Packages](@ref)).

## Scope

**Grid geometries.** An origin (the center of the first voxel), a spacing, a direction matrix for
rotated and reflected grids, and the coordinate space they are given in. Any number of dimensions,
and slices that keep their place in the space of their volume.

**Coordinate spaces as types.** Scanner spaces, the patient spaces LPS and RAS, and spaces that
other packages define. Images and transforms in different spaces are only combined through a
transform that connects them.

**Images with named dimensions.** A `SpatialImage` wraps any array (on the CPU, a GPU or disk) with
its geometry, a name for every dimension, coordinates or labels for the dimensions that are not part
of the grid (time, contrast, coil, echo), and metadata. Grid dimensions may sit anywhere among the
others, as in MPI reconstructions with (contrast, x, y, z, time). Indexing, views, `permutedims`,
`reverse` and broadcasting keep the description.

**Transforms and resampling.** Transforms between spaces are checked when they are composed or
used. Resampling puts an image onto any other grid, through a transform where the spaces differ:
nearest neighbor and multilinear interpolation as one KernelAbstractions kernel for CPUs and GPUs,
B-splines with Interpolations.jl, and resamplers of other packages. Label maps are resampled without
averaging labels.

**Few dependencies and custom components.** Reconstruction packages return these images without
heavy dependencies, and other packages add spaces, resamplers and image types without changing
TomoImages.

**Not in scope.** Reconstruction, registration, segmentation, reading scanner files and displaying
images are tasks of other packages, which exchange their results through these types.

## Where to go next

- [Getting Started](@ref) builds a grid and an image, selects parts of it and resamples it onto a
  rotated grid.
- The tutorials [explore a spatial image](@ref "Explore a Spatial Image") (slices, positions and
  measurements in millimeters, a view along an oblique structure),
  [overlay images from two scanners](@ref "Overlay Images from Two Scanners") through a transform,
  and [analyze a dynamic image](@ref "Analyze a Dynamic Image") with a label map.
- The how-to guides show how to
  [wrap a reconstruction in a SpatialImage](@ref "Wrap a Reconstruction in a SpatialImage"), how to
  [define a space](@ref "Define a Space") and how to [write a resampler](@ref "Write a Resampler").
- [Grids, Spaces and Transforms](@ref) explains the conventions: what the origin, spacing and
  direction mean, what a space is, how dimensions are described, and in which direction transforms
  are stored. [Conventions of Formats and Packages](@ref) compares them with DICOM, NIfTI and other
  packages.
- The reference lists all types and functions ([Spaces and Grids](api/grids.md),
  [Transforms](api/transforms.md), [Images](api/images.md), [Resampling](api/resampling.md)) and the
  [related packages](@ref "Related Packages").
