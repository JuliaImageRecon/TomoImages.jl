"""
    AbstractResampler

Supertype of resamplers: how [`resample`](@ref) interpolates the input at the positions of the output
voxels. A resampler implements

    resample!(r, out, img, map; fillvalue)

(see [`TomoImages.resample!`](@ref)). It may implement [`TomoImages.supports`](@ref) to refuse
inputs before the call, and [`TomoImages.resample_eltype`](@ref) if its output has another element
type than the default. The resampler is the first argument, and its type belongs to the package
that implements it, so that packages add resamplers without defining methods for types of others.

All resamplers of this package read the input over the extent of its voxels: a position up to half
a voxel beyond the outermost voxel centers takes the value at the edge; positions farther out get
`fillvalue`.
"""
abstract type AbstractResampler end

"""
    NearestResampler()

The value of the input voxel nearest to each output position. A position exactly halfway between
two voxel centers takes the voxel with the even index (Julia's `round`). The default for integer
and `Bool` images such as label maps and masks, whose values must not be averaged.
"""
struct NearestResampler <: AbstractResampler end

"""
    LinearResampler()

Multilinear interpolation between the `2^N` input voxels around each output position (bilinear in
2D, trilinear in 3D). Exact at voxel centers, and for images that are linear in the index at
positions between the outermost voxel centers; in the outer half voxel the value of the edge voxel
is taken. The default for images that are not integer or `Bool`.
"""
struct LinearResampler <: AbstractResampler end

"""
    BSplineResampler(degree)

B-spline interpolation of degree 0 to 3 with Interpolations.jl, through the voxel values (an
interpolating spline; degree 2 and 3 prefilter the input, with linear extrapolation as the boundary
condition). Runs on the CPU, for images whose data is an `Array` (not a view or a GPU array).
Throws an `ArgumentError` for other degrees.
"""
struct BSplineResampler <: AbstractResampler
  degree::Int
  function BSplineResampler(degree::Integer)
    0 <= degree <= 3 || throw(ArgumentError("B-spline degree must be 0 to 3, got $degree"))
    new(degree)
  end
end

"""
    resample!(r, out, img, map; fillvalue)

Writes into the array `out` the values of the image `img`, interpolated by the resampler `r`. The
method a custom [`AbstractResampler`](@ref) implements; [`resample`](@ref) allocates `out`,
computes `map` and calls it.

* `out` has the size of the output and lives on the device of `img`; its grid dimensions have the
  size of the target grid, all other dimensions that of `img`.
* For every output index, its indices in the grid dimensions, ordered as the index axes of the
  geometry (`griddims(img)[k]` is the array dimension of index axis `k`), form `j`; `map(j)` (an
  [`TomoImages.IndexMap`](@ref)) is the continuous index of the same position in the grid of `img`,
  in the same order. The indices in all other dimensions are the same in `out` and `img`.
* Positions outside the input get `fillvalue`, which `resample` always passes, converted to the
  element type of `out`.
"""
function resample! end

"""
    default_resampler(img)

The resampler [`resample`](@ref) uses without one: [`NearestResampler`](@ref) for integer and `Bool`
images, [`LinearResampler`](@ref) otherwise. It depends only on the element type; loading a package
never changes it.
"""
default_resampler(img) = eltype(img) <: Union{Integer,Bool} ? NearestResampler() : LinearResampler()

"""
    supports(r, img)

Whether the resampler `r` can resample `img`; [`resample`](@ref) throws an `ArgumentError`
otherwise. Defaults to `true`.
"""
supports(::AbstractResampler, img) = true

"""
    resample([r,] img, target[, t]; fillvalue=zero(resample_eltype(r, eltype(img))))

The image `img` on the geometry `target`: every output voxel takes the value of `img` at its
position, interpolated by the resampler `r` (default: [`TomoImages.default_resampler`](@ref)`(img)`).

* `target`: a [`GridGeometry`](@ref) with as many index and world dimensions as the geometry of
  `img`.
* `t`: without a transform, `target` and the geometry of `img` must be in the same space. With an
  [`AffineTransform`](@ref), positions are mapped from the space of `target` to the space of `img`:
  `source(t)` must be the space of `target` (the fixed image) and `target(t)` that of `img` (the
  moving image). This is the direction registration results are stored in; the output is pulled
  from the input, so every output voxel gets exactly one value.
* `fillvalue`: the value of output voxels whose position lies outside the voxels of `img` (farther
  than half a voxel beyond its outermost voxel centers).

The other dimensions are resampled one by one with the same positions; the description of the
dimensions and the metadata are kept. Throws an `ArgumentError` if the spaces, the number of
dimensions or the resampler (see [`TomoImages.supports`](@ref)) do not fit.

# Example
```julia
g = GridGeometry((64, 64); spacing=(1.0, 1.0), origin=(-31.5, -31.5))
img = TomoImage(rand(Float32, 64, 64), g)
fine = GridGeometry((128, 128); spacing=(0.5, 0.5), origin=(-31.75, -31.75))
resample(img, fine)                          # on a twice finer grid, linear interpolation
resample(NearestResampler(), img, fine; fillvalue=NaN32)
```
"""
resample(img::AbstractTomoImage, target::AbstractGeometry, t=nothing; kwargs...) =
  resample(default_resampler(img), img, target, t; kwargs...)

function resample(r::AbstractResampler, img::AbstractTomoImage, target::GridGeometry, t=nothing;
                  fillvalue=zero(resample_eltype(r, eltype(img))))
  supports(r, img) || throw(ArgumentError("$(nameof(typeof(r))) does not support $(summary(img))"))
  g = geometry(img)
  (ndims(target) == ndims(g) && worlddims(target) == worlddims(g)) ||
    throw(ArgumentError("the target grid has $(ndims(target)) index and $(worlddims(target)) world dimensions, " *
                        "the grid of the image $(ndims(g)) and $(worlddims(g))"))
  if t === nothing
    checkspaces(space(target), space(g), "resample without a transform needs the image in the space of the target grid")
  else
    checkspaces(space(target), source(t), "resample needs a transform from the space of the target grid")
    checkspaces(space(g), TomoImages.target(t), "resample needs a transform into the space of the image")
  end
  data = parent(img)
  Tout = resample_eltype(r, eltype(data))
  imap = index_map(target, g, t; eltype=map_eltype(eltype(data)))
  sd = griddims(img)
  outsize = ntuple(d -> (k = findfirst(==(d), sd); k === nothing ? size(data, d) : size(target, k)), ndims(data))
  out = similar(data, Tout, outsize)
  resample!(r, out, img, imap; fillvalue=convert(Tout, fillvalue))
  TomoImage(out, target, dims(img), meta(img))
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
map_eltype(T) = T <: Union{AbstractFloat,Complex{<:AbstractFloat}} ? promote_type(real(T), Float32) : Float32
