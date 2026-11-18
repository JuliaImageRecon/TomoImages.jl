# Related Packages

TomoImages sits between the packages that reconstruct images and the packages that analyze,
store or show them. This page lists the packages around it and what each of them does.

## Reconstruction

These packages compute images from measurements; their results are what TomoImages describes.

- [CTOperators.jl](https://github.com/JuliaImageRecon/CTOperators.jl): projectors for iterative CT
  reconstruction and FDK, on CPUs and GPUs. Its images are described by an `ImageGeom` of
  ImageGeoms.jl.
- [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl): MRI reconstruction. It
  returns its images as an `AxisArray` with the dimensions (x, y, z, echoes, coils, repetitions).
- [MPIReco.jl](https://github.com/MagneticParticleImaging/MPIReco.jl): reconstruction for magnetic
  particle imaging. It returns an `ImageMeta` around an `AxisArray` with the dimensions (color, x,
  y, z, time), built with ImageUtils.jl.
- [KomaMRI.jl](https://github.com/JuliaHealth/KomaMRI.jl): an MRI simulator. Its phantoms are point
  clouds of spins with tissue properties; it reconstructs through MRIReco.jl.
- [Sinograms.jl](https://github.com/JuliaImageRecon/Sinograms.jl): CT geometries and filtered back
  projection.
- [AbstractImageReconstruction.jl](https://github.com/JuliaImageRecon/AbstractImageReconstruction.jl):
  the framework of reconstruction plans that MPIReco.jl and MRIReco.jl build on.

## Arrays with axes and metadata

These packages attach names, coordinates or metadata to arrays.

- [AxisArrays.jl](https://github.com/JuliaArrays/AxisArrays.jl): one named axis of values per
  dimension; integers are indices, values are selected with `atvalue` and intervals. Each axis is
  independent, so a grid cannot be rotated.
- [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl): named dimensions with lookups,
  selectors (`At`, `Near`, `Between`), stacks of layers and many integrations. Rotated grids use a
  `Transformed` lookup.
- [Rasters.jl](https://github.com/rafaqz/Rasters.jl): builds on DimensionalData.jl and adds
  coordinate reference systems, reprojection and resampling for geospatial data. It does for maps
  roughly what TomoImages does for tomographic images.
- [ImageCore.jl](https://github.com/JuliaImages/ImageCore.jl) and the JuliaImages packages around
  it: traits such as `pixelspacing` and `spacedirections` for images, colors, and image
  processing. There is no trait for the origin of an image.
- [ImageGeoms.jl](https://github.com/JuliaImageRecon/ImageGeoms.jl): the image grid of iterative
  reconstruction (dimensions, pixel sizes, offsets in pixels from the center, a support mask), with
  the conventions of the Michigan Image Reconstruction Toolbox.
- [MedImages.jl](https://github.com/JuliaHealth/MedImages.jl): 3D medical images with origin,
  spacing and direction, and patient metadata.

TomoImages differs from AxisArrays.jl and DimensionalData.jl in what it puts first: one geometry for
all grid dimensions together (origin, spacing and a direction matrix, so rotated grids and slices
are exact), the coordinate space as a type that is checked, transforms between spaces, and
resampling between grids. Its non-grid dimensions are deliberately simple; for analysis along them,
the packages above offer far more.

## Algorithms and formats

- [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl): the B-spline interpolation
  behind [`BSplineResampler`](@ref).
- [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl): the kernels of
  [`NearestResampler`](@ref) and [`LinearResampler`](@ref), for CPUs and GPUs.
- [ImageTransformations.jl](https://github.com/JuliaImages/ImageTransformations.jl): warping images
  in index space, with
  [CoordinateTransformations.jl](https://github.com/JuliaGeometry/CoordinateTransformations.jl) for
  the transforms.
- [NIfTI.jl](https://github.com/JuliaNeuroscience/NIfTI.jl) and
  [DICOM.jl](https://github.com/JuliaHealth/DICOM.jl): reading and writing medical image files.

[Conventions of Formats and Packages](@ref) shows how the grids of DICOM, NIfTI, ImageGeoms.jl and
the reconstruction packages translate into grid geometries.
