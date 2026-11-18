# Selections that keep the geometry: integers drop a dimension, ranges and colons keep it.
const SliceIndex = Union{Integer,AbstractRange{<:Integer},Colon}

Base.getindex(img::SpatialImage{T,N}, I::Vararg{SliceIndex,N}) where {T,N} = select(getindex, img, I)
Base.view(img::SpatialImage{T,N}, I::Vararg{SliceIndex,N}) where {T,N} = select(view, img, I)

# Selection by name, with indices: img[x=10:20, time=3]. Unnamed dimensions are kept whole.
function Base.getindex(img::SpatialImage; kwargs...)
  isempty(kwargs) && return img
  img[named_indices(img, kwargs)...]
end

function named_indices(img, kwargs)
  names = keys(dims(img))
  for k in keys(kwargs)
    k in names || throw(ArgumentError("the image has no dimension $k; its dimensions are $(join(names, ", "))"))
    kwargs[k] isa SliceIndex ||
      throw(ArgumentError("dimension $k: indices are integers, ranges or colons, got $(repr(kwargs[k])); " *
                          "select by value with selectvalue"))
  end
  ntuple(d -> get(kwargs, names[d], Colon()), ndims(img))
end

"""
    selectvalue(img; name=value, ...)

The part of `img` at the given values of its named non-grid dimensions: the nearest coordinate on a
numeric axis, the equal label on a [`CategoricalAxis`](@ref) (see [`TomoImages.indexof`](@ref)).
Each selected dimension is dropped. Indexing with integers and ranges, `img[time=3]`, always means
indices.

# Example
```julia
selectvalue(img; contrast="tracer B", time=1.5)
```
"""
function selectvalue(img::AbstractSpatialImage; kwargs...)
  D = dims(img)
  names = keys(D)
  for k in keys(kwargs)
    k in names || throw(ArgumentError("the image has no dimension $k; its dimensions are $(join(names, ", "))"))
    D[k] isa GridAxis &&
      throw(ArgumentError("dimension $k is a grid axis, which has no values of its own; use index_to_world " *
                          "and world_to_index, or resample"))
  end
  img[ntuple(d -> haskey(kwargs, names[d]) ? indexof(D[d], kwargs[names[d]]) : Colon(), ndims(img))...]
end

function select(f, img::SpatialImage, I::Tuple)
  N = ndims(img)
  J = ntuple(d -> I[d] isa Colon ? axes(img, d) : I[d], N)
  data = f(img.data, J...)
  sd = spatialdims(img)
  all(k -> J[sd[k]] isa Integer, eachindex(sd)) && return data
  geom = subgrid(img.geom, map(d -> J[d], sd))
  # grid axes that remain are numbered in the order of the index axes of the original geometry
  keptaxes = [k for k in eachindex(sd) if !(J[sd[k]] isa Integer)]
  D = img.dims
  kept = [d for d in 1:N if !(J[d] isa Integer)]
  newdims = NamedTuple{Tuple(keys(D)[d] for d in kept)}(Tuple(
    D[d] isa GridAxis ? GridAxis(findfirst(==(D[d].index), keptaxes)) : D[d][J[d]] for d in kept))
  SpatialImage(data, geom, newdims, img.meta)
end

function Base.permutedims(img::SpatialImage, perm)
  p = Tuple(perm)
  newdims = NamedTuple{map(d -> keys(img.dims)[d], p)}(map(d -> img.dims[d], p))
  SpatialImage(permutedims(img.data, p), img.geom, newdims, img.meta)
end

function Base.reverse(img::SpatialImage; dims=:)
  ds = dims isa Colon ? (1:ndims(img)) : dims
  img[ntuple(d -> d in ds ? (size(img, d):-1:1) : Colon(), ndims(img))...]
end
