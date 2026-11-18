"""
    AbstractSpatialImage{T,N}

Supertype of arrays that know where their voxels are. An image implements the `AbstractArray`
interface and the accessors `parent`, [`geometry`](@ref), [`TomoImages.dims`](@ref) and
[`meta`](@ref); functions of this package use only these, never the fields of an image type.
"""
abstract type AbstractSpatialImage{T,N} <: AbstractArray{T,N} end

"""
    SpatialImage(data, geometry; dims, meta=(;))

An array `data` whose grid dimensions lie on `geometry`, with a named description of every
dimension.

* `dims`: a `NamedTuple` with one entry per dimension of `data`, in order. A [`GridAxis`](@ref)`(k)`
  marks index axis `k` of the geometry; every grid axis appears once. Any other dimension is
  described by a vector with one value per index: a [`NumericAxis`](@ref) (ordered coordinates with
  a unit), a [`CategoricalAxis`](@ref) (labels), or any `AbstractVector`. The names are part of the
  type. Default: the grid axes `x1, x2, ...` first, then `dim_d = 1:n` for further dimensions.
* `meta`: free metadata; geometry code never reads it.

Indexing, `view`, `permutedims`, `reverse` and broadcasting keep the geometry, the axes and the
metadata wherever the result still has grid axes. `parent(img)` is the array, on whatever device it
lives.

# Example
```julia
g = GridGeometry((40, 40, 20); spacing=(1.0, 1.0, 2.0))
img = SpatialImage(rand(Float32, 2, 40, 40, 20, 100), g;
                   dims=(contrast=CategoricalAxis(["tracer A", "tracer B"]),
                         x=GridAxis(1), y=GridAxis(2), z=GridAxis(3),
                         time=NumericAxis(range(0, step=0.02, length=100); unit="s")))
img[contrast=2, time=51]                    # by index: a 40x40x20 SpatialImage on g
selectvalue(img; contrast="tracer B", time=1.0)   # the same, by value
```
"""
struct SpatialImage{T,N,A<:AbstractArray{T,N},G<:AbstractGeometry,D<:NamedTuple,M} <: AbstractSpatialImage{T,N}
  data::A
  geom::G
  dims::D
  meta::M
end

function SpatialImage(data::AbstractArray, geom::AbstractGeometry; dims::NamedTuple=default_dims(Val(ndims(geom)), ndims(data)),
                      meta=NamedTuple())
  N, Ns = ndims(data), ndims(geom)
  length(dims) == N || throw(ArgumentError("dims needs one entry per dimension of the data ($N), got $(length(dims))"))
  # fill the default index vectors of non-grid dimensions
  dims = NamedTuple{keys(dims)}(ntuple(d -> dims[d] === nothing ? Base.OneTo(size(data, d)) : dims[d], N))
  grid = [d for d in 1:N if dims[d] isa GridAxis]
  sort([dims[d].index for d in grid]) == 1:Ns ||
    throw(ArgumentError("every index axis 1 to $Ns of the geometry needs exactly one GridAxis in dims, got " *
                        "$(Tuple(dims[d].index for d in grid))"))
  for d in 1:N
    name, a = keys(dims)[d], dims[d]
    if a isa GridAxis
      size(data, d) == size(geom, a.index) ||
        throw(DimensionMismatch("dimension $d ($name) of the data has size $(size(data, d)), index axis " *
                                "$(a.index) of the geometry $(size(geom, a.index))"))
    elseif a isa AbstractVector
      length(a) == size(data, d) ||
        throw(DimensionMismatch("axis $name has length $(length(a)), dimension $d of the data $(size(data, d))"))
    else
      throw(ArgumentError("axis $name must be a GridAxis or an AbstractVector, got $(typeof(a))"))
    end
  end
  SpatialImage(data, geom, dims, meta)
end

Base.parent(img::SpatialImage) = img.data

"""
    geometry(img)

The geometry of the grid dimensions of `img`.
"""
geometry(img::SpatialImage) = img.geom

"""
    dims(img)

The description of every dimension of `img`, in order, as a `NamedTuple`: a [`GridAxis`](@ref) for
the grid dimensions and a vector of values for the others. Public, not exported: DimensionalData.jl
exports a function of the same name.
"""
dims(img::SpatialImage) = img.dims

"""
    meta(img)

The metadata of `img`.
"""
meta(img::SpatialImage) = img.meta

"""
    spatialdims(img)

The array dimensions of the index axes `1, 2, ...` of the geometry of `img`: `spatialdims(img)[k]`
is the dimension marked `GridAxis(k)`.
"""
function spatialdims(img::AbstractSpatialImage)
  D = dims(img)
  ntuple(k -> findfirst(a -> a isa GridAxis && a.index == k, Tuple(D)), ndims(geometry(img)))
end

space(img::AbstractSpatialImage) = space(geometry(img))

# the same image description around other data of the same shape
rewrap(img::SpatialImage, data) = SpatialImage(data, img.geom, img.dims, img.meta)

Base.size(img::SpatialImage) = size(img.data)
Base.IndexStyle(::Type{<:SpatialImage{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.@propagate_inbounds Base.getindex(img::SpatialImage{T,N}, I::Vararg{Int,N}) where {T,N} = img.data[I...]
Base.@propagate_inbounds Base.getindex(img::SpatialImage, i::Int) = img.data[i]
Base.@propagate_inbounds Base.setindex!(img::SpatialImage{T,N}, v, I::Vararg{Int,N}) where {T,N} = (img.data[I...] = v; img)
Base.@propagate_inbounds Base.setindex!(img::SpatialImage, v, i::Int) = (img.data[i] = v; img)

Base.similar(img::SpatialImage, ::Type{S}, dims::Dims) where {S} =
  dims == size(img) ? rewrap(img, similar(img.data, S)) : similar(img.data, S, dims)
Base.copy(img::SpatialImage) = rewrap(img, copy(img.data))

Adapt.adapt_structure(to, img::SpatialImage) = rewrap(img, adapt(to, img.data))

function Base.summary(io::IO, img::SpatialImage)
  print(io, Base.dims2string(size(img)), " SpatialImage{", eltype(img), "} in ", spacename(space(img)),
        " with dims ", join(keys(img.dims), ", "))
end
