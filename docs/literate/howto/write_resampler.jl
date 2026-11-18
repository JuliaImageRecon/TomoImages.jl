# # Write a Resampler
#
# [`resample`](@ref) takes the resampler as its first argument. A resampler is a type of its own, so
# a package adds one by defining a type and a method for it, without changing TomoImages. This guide
# writes a resampler for masks that keeps thin structures, which nearest-neighbor resampling loses
# on a coarser grid.

using TomoImages

# ## The problem
#
# A mask of a vessel one pixel wide, on a grid of 1 mm:

g = GridGeometry((32, 32); spacing=(1.0, 1.0), origin=(0.0, 0.0))
vessel = SpatialImage([i == 15 for i in 1:32, j in 1:32], g);

# On a grid of 2 mm, every coarse pixel lies halfway between two columns of fine pixels.
# [`NearestResampler`](@ref) reads one of them, rounding halfway positions to the even index, and
# here never reads the column of the vessel: the coarse mask is empty.

coarse = GridGeometry((16, 16); spacing=(2.0, 2.0), origin=(0.5, 0.5))
count(parent(resample(NearestResampler(), vessel, coarse)))

# ## The resampler
#
# `AnyNeighbor` marks an output pixel if any of the input pixels around its position is marked.
# A resampler implements [`resample!`](@ref): it fills the output array `out` with the values of the
# image at the continuous input indices `imap(j)` of the output indices `j`. This one supports images
# without extra axes on the CPU, which [`supports`](@ref) states, and returns `Bool`, which
# [`TomoImages.resample_eltype`](@ref) states:

struct AnyNeighbor <: AbstractResampler end

TomoImages.supports(::AnyNeighbor, img) = eltype(img) == Bool && parent(img) isa Array && isempty(extra(img))
TomoImages.resample_eltype(::AnyNeighbor, ::Type) = Bool

function TomoImages.resample!(::AnyNeighbor, out, img, imap; fill=false)
  data = parent(img)
  n = size(data)
  for j in CartesianIndices(out)
    p = imap(Tuple(j))
    if all(k -> 0.5 <= p[k] <= n[k] + 0.5, eachindex(p))
      lo = ntuple(k -> clamp(floor(Int, p[k]), 1, n[k]), length(n))
      hi = ntuple(k -> clamp(floor(Int, p[k]) + 1, 1, n[k]), length(n))
      out[j] = any(data[CartesianIndex(lo):CartesianIndex(hi)])
    else
      out[j] = fill
    end
  end
  out
end

# The same call as with the built-in resamplers; every coarse pixel along the vessel is marked:

mask = resample(AnyNeighbor(), vessel, coarse)
count(parent(mask)), eltype(mask)

# For an image it does not support, `resample` refuses before the call:

try #hide
  resample(AnyNeighbor(), SpatialImage(rand(32, 32), g), coarse)
catch e #hide
  print(sprint(showerror, e)) #hide
end #hide

# The index map `imap` already contains the two geometries and the transform, if any, so a
# resampler never handles spaces or directions. On a GPU, a resampler runs its loop as a kernel; the
# map is `isbits` and is passed as an argument (see [`TomoImages.IndexMap`](@ref)).
