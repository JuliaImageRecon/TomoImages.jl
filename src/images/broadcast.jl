# Broadcasting keeps the geometry of the images among the arguments, which must all be in the same
# space and on the same grid. The computation itself runs on the parents, so that GPU arrays use
# their own broadcast kernels instead of indexing the wrapper element by element.
struct TomoStyle <: Broadcast.AbstractArrayStyle{Any} end
TomoStyle(::Val) = TomoStyle()

Base.BroadcastStyle(::Type{<:TomoImage}) = TomoStyle()
Base.BroadcastStyle(s::TomoStyle, ::Broadcast.AbstractArrayStyle) = s
Base.BroadcastStyle(s::TomoStyle, ::Broadcast.DefaultArrayStyle) = s
Base.BroadcastStyle(s::TomoStyle, ::TomoStyle) = s

firstimage(bc::Broadcast.Broadcasted) = firstimage(bc.args)
firstimage(args::Tuple) = firstimage(firstimage(args[1]), Base.tail(args))
firstimage(img::TomoImage, rest) = img
firstimage(::Any, rest) = firstimage(rest)
firstimage(x) = x
firstimage(::Tuple{}) = nothing

unwrap(bc::Broadcast.Broadcasted) = Broadcast.broadcasted(bc.f, map(unwrap, bc.args)...)
unwrap(img::TomoImage) = parent(img)
unwrap(x) = x

images(bc::Broadcast.Broadcasted) = mapreduce(images, vcat, bc.args; init=AbstractTomoImage[])
images(img::AbstractTomoImage) = AbstractTomoImage[img]
images(x) = AbstractTomoImage[]

function checkgrids(bc, dest=nothing)
  imgs = dest === nothing ? images(bc) : vcat(images(dest), images(bc))
  for img in imgs[2:end]
    checkspaces(space(imgs[1]), space(img), "broadcasting images needs them in one space")
    isapprox(geometry(imgs[1]), geometry(img)) ||
      throw(ArgumentError("broadcasting images needs them on the same grid; resample one onto the grid of the other"))
  end
  nothing
end

function Base.similar(bc::Broadcast.Broadcasted{TomoStyle}, ::Type{T}) where {T}
  checkgrids(bc)
  img = firstimage(bc)
  data = similar(parent(img), T, axes(bc))
  axes(bc) == axes(img) ? rewrap(img, data) : data
end

function Base.copyto!(dest::TomoImage, bc::Broadcast.Broadcasted{TomoStyle})
  checkgrids(bc, dest)
  copyto!(parent(dest), Broadcast.instantiate(unwrap(bc)))
  dest
end
