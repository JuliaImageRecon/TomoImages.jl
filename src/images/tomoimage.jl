"""
    AbstractTomoImage{T,N}

Supertype of arrays that know where their voxels are. An image implements the `AbstractArray`
interface and the accessors `parent`, [`geometry`](@ref), [`TomoImages.dims`](@ref) and
[`meta`](@ref); functions of this package use only these, never the fields of an image type.
"""
abstract type AbstractTomoImage{T,N} <: AbstractArray{T,N} end

"""
    TomoImage(data, geom; dims=(x1=GridAxis(1), ..., dim_d=1:n, ...), meta=(;))

A tomographic image: an array `data` whose grid dimensions lie on the geometry `geom`, with a named
description of every dimension.

* `dims`: a `NamedTuple` with one entry per dimension of `data`, in order. A [`GridAxis`](@ref)`(k)`
  marks index axis `k` of the geometry; every grid axis appears once. Any other dimension is
  described by a vector with one value per index: a [`NumericAxis`](@ref) (ordered coordinates with
  a unit), a [`CategoricalAxis`](@ref) (labels), or any `AbstractVector`. The names are part of the
  type. Default: the grid axes `x1, x2, ...` first, then `dim_d = 1:n` for further dimensions.
* `meta`: free metadata; geometry code never reads it.

Indexing and `view` with integers, ranges and colons (also by name, `img[time=3]`), `permutedims`,
`reverse` and broadcasting keep the geometry, the dimensions and the metadata wherever the result
still has grid dimensions; other indices (vectors, logical masks) give a plain array. Broadcasting
several images requires them in one space and on the same grid. `parent(img)` is the array, on
whatever device it lives.

Throws an `ArgumentError` or a `DimensionMismatch` if `dims` does not fit the data and the geometry.

# Example
```julia
g = GridGeometry((40, 40, 20); spacing=(1.0, 1.0, 2.0))
img = TomoImage(rand(Float32, 2, 40, 40, 20, 100), g;
                   dims=(contrast=CategoricalAxis(["tracer A", "tracer B"]),
                         x=GridAxis(1), y=GridAxis(2), z=GridAxis(3),
                         time=NumericAxis(range(0, step=0.02, length=100); unit="s")))
img[contrast=2, time=51]                    # by index: a 40x40x20 TomoImage on g
selectvalue(img; contrast="tracer B", time=1.0)   # the same, by value
```
"""
struct TomoImage{T,N,A<:AbstractArray{T,N},G<:AbstractGeometry,D<:NamedTuple,M} <: AbstractTomoImage{T,N}
  data::A
  geom::G
  dims::D
  meta::M
end

function TomoImage(data::AbstractArray, geom::AbstractGeometry;
                   dims::NamedTuple=default_dims(Val(ndims(geom)), ndims(data)),
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
  TomoImage(data, geom, dims, meta)
end

Base.parent(img::TomoImage) = img.data

"""
    geometry(img)

The geometry of the grid dimensions of `img`.
"""
geometry(img::TomoImage) = img.geom

"""
    dims(img)

The description of every dimension of `img`, in order, as a `NamedTuple`: a [`GridAxis`](@ref) for
the grid dimensions and a vector of values for the others. Public, not exported: DimensionalData.jl
exports a function of the same name.
"""
dims(img::TomoImage) = img.dims

"""
    meta(img)

The metadata of `img`.
"""
meta(img::TomoImage) = img.meta

"""
    griddims(img)

The array dimensions of the index axes `1, 2, ...` of the geometry of `img`: `griddims(img)[k]`
is the dimension marked `GridAxis(k)`.
"""
function griddims(img::AbstractTomoImage)
  D = dims(img)
  ntuple(k -> findfirst(a -> a isa GridAxis && a.index == k, Tuple(D)), ndims(geometry(img)))
end

space(img::AbstractTomoImage) = space(geometry(img))

# the same image description around other data of the same shape
rewrap(img::TomoImage, data) = TomoImage(data, img.geom, img.dims, img.meta)

Base.size(img::TomoImage) = size(img.data)
Base.IndexStyle(::Type{<:TomoImage{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.@propagate_inbounds Base.getindex(img::TomoImage{T,N}, I::Vararg{Int,N}) where {T,N} = img.data[I...]
Base.@propagate_inbounds Base.getindex(img::TomoImage, i::Int) = img.data[i]
Base.@propagate_inbounds Base.setindex!(img::TomoImage{T,N}, v, I::Vararg{Int,N}) where {T,N} =
  (img.data[I...] = v; img)
Base.@propagate_inbounds Base.setindex!(img::TomoImage, v, i::Int) = (img.data[i] = v; img)

Base.similar(img::TomoImage, ::Type{S}, dims::Dims) where {S} =
  dims == size(img) ? rewrap(img, similar(img.data, S)) : similar(img.data, S, dims)
Base.copy(img::TomoImage) = rewrap(img, copy(img.data))

Adapt.adapt_structure(to, img::TomoImage) = rewrap(img, adapt(to, img.data))

function Base.summary(io::IO, img::TomoImage)
  print(io, Base.dims2string(size(img)), " TomoImage{", eltype(img), "} in ", spacename(space(img)),
        " with dims ", join(keys(img.dims), ", "))
end
