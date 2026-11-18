# Selections that keep the geometry: integers drop a dimension, ranges and colons keep it.
const SliceIndex = Union{Integer,AbstractRange{<:Integer},Colon}

Base.getindex(img::SpatialImage{T,N}, I::Vararg{SliceIndex,N}) where {T,N} = select(getindex, img, I)
Base.view(img::SpatialImage{T,N}, I::Vararg{SliceIndex,N}) where {T,N} = select(view, img, I)

# Selection by value on the extra axes: img[time=1.5, contrast="tracer A"]
function Base.getindex(img::SpatialImage; kwargs...)
  isempty(kwargs) && return img
  names = keys(img.extra)
  for k in keys(kwargs)
    k in names || throw(ArgumentError("the image has no extra axis $k; its extra axes are $(join(names, ", "))"))
  end
  others = Tuple(d for d in 1:ndims(img) if !(d in img.spatialdims))
  I = ntuple(ndims(img)) do d
    k = findfirst(==(d), others)
    k === nothing && return Colon()
    name = names[k]
    haskey(kwargs, name) ? indexof(img.extra[name], kwargs[name]) : Colon()
  end
  img[I...]
end

function select(f, img::SpatialImage, I::Tuple)
  N = ndims(img)
  J = ntuple(d -> I[d] isa Colon ? axes(img, d) : I[d], N)
  data = f(img.data, J...)
  sd = img.spatialdims
  all(d -> J[d] isa Integer, sd) && return data
  geom = subgrid(img.geom, map(d -> J[d], sd))
  kept = [d for d in 1:N if !(J[d] isa Integer)]
  newsd = Tuple(findfirst(==(d), kept) for d in sd if !(J[d] isa Integer))
  others = [d for d in 1:N if !(d in sd)]
  names = keys(img.extra)
  pairs_kept = [names[k] => img.extra[k][J[d]] for (k, d) in enumerate(others) if !(J[d] isa Integer)]
  SpatialImage(data, geom, newsd, NamedTuple(pairs_kept), img.meta)
end
