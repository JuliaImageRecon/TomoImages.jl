"""
    AbstractSpatialImage{T,N}

Supertype of arrays that know where their voxels are. An image implements the `AbstractArray`
interface and the accessors `parent`, [`geometry`](@ref), [`spatialdims`](@ref), [`extra`](@ref) and
[`meta`](@ref); functions of this package use only these, never the fields of an image type.
"""
abstract type AbstractSpatialImage{T,N} <: AbstractArray{T,N} end

"""
    SpatialImage(data, geometry; spatialdims, extra=(;), meta=(;))

An array `data` whose spatial dimensions lie on `geometry`, with named axes for the other
dimensions.

* `spatialdims`: the dimensions of `data` that are spatial, in the order of the index axes of the
  geometry, increasing (default: the first `ndims(geometry)` dimensions). MPI data with dimensions
  (contrast, x, y, z, time) has `spatialdims = (2, 3, 4)`.
* `extra`: a `NamedTuple` with one [`AbstractAxis`](@ref) for every other dimension, in order. The
  names are part of the type.
* `meta`: free metadata; geometry code never reads it.

Indexing, `view` and broadcasting keep the geometry, the axes and the metadata wherever the result
still has spatial dimensions (see [`subgrid`](@ref)). `parent(img)` is the array, on whatever device
it lives.

# Example
```julia
g = GridGeometry((40, 40, 20); spacing=(1.0, 1.0, 2.0))
img = SpatialImage(rand(Float32, 2, 40, 40, 20, 100), g; spatialdims=(2, 3, 4),
                   extra=(contrast=CategoricalAxis(["tracer A", "tracer B"]),
                          time=NumericAxis(range(0, step=0.02, length=100); unit="s")))
img[contrast="tracer B", time=1.0]         # a 40x40x20 SpatialImage on g
```
"""
struct SpatialImage{T,N,A<:AbstractArray{T,N},G<:AbstractGeometry,P,X<:NamedTuple,M} <: AbstractSpatialImage{T,N}
  data::A
  geom::G
  spatialdims::P
  extra::X
  meta::M
end

function SpatialImage(data::AbstractArray, geom::AbstractGeometry; spatialdims=ntuple(identity, ndims(geom)),
                      extra::NamedTuple=NamedTuple(), meta=NamedTuple())
  N, Ns = ndims(data), ndims(geom)
  sd = Tuple(Int.(spatialdims))
  length(sd) == Ns || throw(ArgumentError("the geometry has $Ns dimensions, spatialdims names $(length(sd))"))
  (issorted(sd) && allunique(sd) && all(d -> 1 <= d <= N, sd)) ||
    throw(ArgumentError("spatialdims must be increasing dimensions of the data (1 to $N), got $sd"))
  all(k -> size(data, sd[k]) == size(geom, k), 1:Ns) ||
    throw(DimensionMismatch("the spatial dimensions $sd of the data have size $(map(d -> size(data, d), sd)), " *
                            "the geometry $(size(geom))"))
  others = Tuple(d for d in 1:N if !(d in sd))
  length(extra) == length(others) ||
    throw(ArgumentError("the data has $(length(others)) non-spatial dimensions, extra names $(length(extra))"))
  for (k, (name, axis)) in enumerate(pairs(extra))
    axis isa AbstractAxis || throw(ArgumentError("extra axis $name must be an AbstractAxis, got $(typeof(axis))"))
    length(axis) == size(data, others[k]) ||
      throw(DimensionMismatch("extra axis $name has length $(length(axis)), dimension $(others[k]) of the data " *
                              "$(size(data, others[k]))"))
  end
  SpatialImage(data, geom, sd, extra, meta)
end

Base.parent(img::SpatialImage) = img.data

"""
    geometry(img)

The geometry of the spatial dimensions of `img`.
"""
geometry(img::SpatialImage) = img.geom

"""
    spatialdims(img)

The dimensions of `img` that are spatial, in the order of the index axes of its geometry.
"""
spatialdims(img::SpatialImage) = img.spatialdims

"""
    extra(img)

The named axes of the non-spatial dimensions of `img`, in order, as a `NamedTuple`.
"""
extra(img::SpatialImage) = img.extra

"""
    meta(img)

The metadata of `img`.
"""
meta(img::SpatialImage) = img.meta

space(img::AbstractSpatialImage) = space(geometry(img))

# the same image description around other data of the same shape
rewrap(img::SpatialImage, data) = SpatialImage(data, img.geom, img.spatialdims, img.extra, img.meta)

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
  print(io, Base.dims2string(size(img)), " SpatialImage{", eltype(img), "} in ", spacename(space(img)))
  isempty(img.extra) || print(io, ", extra axes ", join(keys(img.extra), ", "))
end
