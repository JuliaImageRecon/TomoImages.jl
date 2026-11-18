"""
    AbstractResampler

Supertype of resamplers: how [`resample`](@ref) interpolates the input at the positions of the output
voxels. A resampler implements

    resample!(r, out, img, map)

which writes into `out` (an array of the output size, on the device of `img`) the values of `img` at
the continuous input indices `map(j)` of the output spatial indices `j` (an [`IndexMap`](@ref)),
with every non-spatial index unchanged. It may implement [`supports`](@ref) to refuse inputs before
the call, and [`TomoImages.resample_eltype`](@ref) if its output has another element type than the
default. The resampler is the first argument, and its type belongs to the package that implements
it, so that packages add resamplers without defining methods for types of others.

All resamplers of this package read the input over the extent of its voxels: a position up to half
a voxel beyond the outermost voxel centers takes the value at the edge; positions farther out get
`fill`.
"""
abstract type AbstractResampler end

"""
    NearestResampler()

The value of the input voxel nearest to each output position. The default for integer images such
as label maps, whose values are names that must not be averaged.
"""
struct NearestResampler <: AbstractResampler end

"""
    LinearResampler()

Multilinear interpolation between the `2^N` input voxels around each output position (bilinear in
2D, trilinear in 3D). Exact at voxel centers and for images that are linear in the index. The
default for real and complex images.
"""
struct LinearResampler <: AbstractResampler end

"""
    BSplineResampler(degree)

B-spline interpolation of degree 0 to 3 with Interpolations.jl, through the voxel values (an
interpolating spline; degree 2 and 3 prefilter the input). Runs on the CPU, for images whose data
is an `Array`.
"""
struct BSplineResampler <: AbstractResampler
  degree::Int
  function BSplineResampler(degree::Integer)
    0 <= degree <= 3 || throw(ArgumentError("B-spline degree must be 0 to 3, got $degree"))
    new(degree)
  end
end

"""
    resample!(r, out, img, map; fill)

Writes into the array `out` the values of the image `img` at the continuous input indices `map(j)`
of the output spatial indices `j`, interpolated by the resampler `r`; positions outside the input
get `fill`. `out` has the size of the output and lives on the device of `img`; `map` is an
[`TomoImages.IndexMap`](@ref). The method a custom [`AbstractResampler`](@ref) implements;
[`resample`](@ref) allocates `out`, computes `map` and calls it.
"""
function resample! end

"""
    default_resampler(img)

The resampler [`resample`](@ref) uses without one: [`NearestResampler`](@ref) for integer and `Bool` images,
[`LinearResampler`](@ref) otherwise. Loading a package never changes it.
"""
default_resampler(img) = eltype(img) <: Union{Integer,Bool} ? NearestResampler() : LinearResampler()

"""
    supports(r, img)

Whether the resampler `r` can resample `img`; [`resample`](@ref) refuses otherwise. Defaults to `true`.
"""
supports(::AbstractResampler, img) = true

"""
    resample([r,] img, target[, t]; fill=zero)

The image `img` on the geometry `target`: every output voxel takes the value of `img` at its
position, interpolated by the resampler `r` (default: [`default_resampler`](@ref)`(img)`).

Without a transform, `target` and the geometry of `img` must be in the same space. With an
[`AffineTransform`](@ref) `t`, positions are mapped from the space of `target` to the space of `img`:
`source(t)` must be the space of `target` (the fixed image) and `target(t)` that of `img` (the moving
image). This is the direction registration results are stored in; the output is pulled from the
input, so every output voxel gets exactly one value.

Non-spatial dimensions are resampled one by one with the same positions; the axes and the metadata
are kept. Positions outside the input get `fill` (default zero).

# Example
```julia
fine = GridGeometry((128, 128); spacing=(0.5, 0.5), origin=(-31.75, -31.75))
resample(img, fine)                          # img on a twice finer grid, linear interpolation
resample(NearestResampler(), labels, fine)
```
"""
resample(img::AbstractSpatialImage, target::AbstractGeometry, t=nothing; kwargs...) =
  resample(default_resampler(img), img, target, t; kwargs...)

function resample(r::AbstractResampler, img::AbstractSpatialImage, target::GridGeometry, t=nothing;
                  fill=nothing)
  supports(r, img) || throw(ArgumentError("$(typeof(r)) does not support $(summary(img))"))
  g = geometry(img)
  ndims(target) == ndims(g) ||
    throw(ArgumentError("the target grid has $(ndims(target)) dimensions, the image $(ndims(g)) spatial ones"))
  if t === nothing
    checkspaces(space(target), space(g), "resample without a transform needs the image in the space of the target grid")
  else
    checkspaces(space(target), source(t), "resample needs a transform from the space of the target grid")
    checkspaces(space(g), TomoImages.target(t), "resample needs a transform into the space of the image")
  end
  data = parent(img)
  Tout = resample_eltype(r, eltype(data))
  imap = index_map(target, g, t; eltype=map_eltype(eltype(data)))
  sd = spatialdims(img)
  outsize = ntuple(d -> (k = findfirst(==(d), sd); k === nothing ? size(data, d) : size(target, k)), ndims(data))
  out = similar(data, Tout, outsize)
  resample!(r, out, img, imap; fill=fill === nothing ? zero(Tout) : convert(Tout, fill))
  SpatialImage(out, target, dims(img), meta(img))
end

"""
    resample_eltype(r, T)

The element type of the output of the resampler `r` for input elements of type `T`. Interpolation
averages, so by default integers and `Bool` become `Float32` and other types stay;
[`NearestResampler`](@ref) keeps every type. A resampler that returns another type implements a
method.
"""
resample_eltype(::AbstractResampler, T) = T <: Union{Integer,Bool} ? Float32 : T
resample_eltype(::NearestResampler, T) = T

# The precision of the index map: that of the data, at least Float32 (GPUs without Float64 run the
# kernels in Float32).
map_eltype(T) = T <: Union{AbstractFloat,Complex{<:AbstractFloat}} ? real(T) : Float32
