# Broadcasting keeps the geometry of the first image among the arguments. The computation itself
# runs on the parents, so that GPU arrays use their own broadcast kernels instead of indexing the
# wrapper element by element.
struct SpatialStyle <: Broadcast.AbstractArrayStyle{Any} end
SpatialStyle(::Val) = SpatialStyle()

Base.BroadcastStyle(::Type{<:SpatialImage}) = SpatialStyle()
Base.BroadcastStyle(s::SpatialStyle, ::Broadcast.AbstractArrayStyle) = s
Base.BroadcastStyle(s::SpatialStyle, ::Broadcast.DefaultArrayStyle) = s
Base.BroadcastStyle(s::SpatialStyle, ::SpatialStyle) = s

firstimage(bc::Broadcast.Broadcasted) = firstimage(bc.args)
firstimage(args::Tuple) = firstimage(firstimage(args[1]), Base.tail(args))
firstimage(img::SpatialImage, rest) = img
firstimage(::Any, rest) = firstimage(rest)
firstimage(x) = x
firstimage(::Tuple{}) = nothing

unwrap(bc::Broadcast.Broadcasted) = Broadcast.broadcasted(bc.f, map(unwrap, bc.args)...)
unwrap(img::SpatialImage) = parent(img)
unwrap(x) = x

function Base.similar(bc::Broadcast.Broadcasted{SpatialStyle}, ::Type{T}) where {T}
  img = firstimage(bc)
  data = similar(parent(img), T, axes(bc))
  axes(bc) == axes(img) ? rewrap(img, data) : data
end

function Base.copyto!(dest::SpatialImage, bc::Broadcast.Broadcasted{SpatialStyle})
  copyto!(parent(dest), Broadcast.instantiate(unwrap(bc)))
  dest
end
