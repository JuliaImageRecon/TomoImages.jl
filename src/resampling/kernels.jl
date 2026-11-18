# Nearest-neighbor and multilinear resampling as one KernelAbstractions kernel for CPUs and GPUs: one work
# item per output element. The grid dimensions are given by `sd` (their positions in the data)
# and `pos` (for every dimension of the data its index among the grid ones, 0 for the others),
# both as type parameters, so that index tuples are built at compile time.

@inline grid_index(I::NTuple, ::Val{SD}) where {SD} = ntuple(k -> I[SD[k]], Val(length(SD)))

@inline with_grid(I::NTuple{N,Int}, ::Val{POS}, s) where {N,POS} = ntuple(d -> POS[d] == 0 ? I[d] : s[POS[d]], Val(N))

# inside the extent of the input voxels: up to half a voxel beyond the outermost centers
@inline inside(p::SVector{Ns,T}, n::NTuple{Ns,Int}) where {Ns,T} =
  all(ntuple(k -> T(0.5) <= p[k] <= n[k] + T(0.5), Val(Ns)))

@inline function sample(::NearestResampler, data, I, pos, n::NTuple{Ns,Int}, p, fillvalue) where {Ns}
  inside(p, n) || return fillvalue
  idx = ntuple(k -> clamp(unsafe_trunc(Int, round(p[k])), 1, n[k]), Val(Ns))
  @inbounds data[with_grid(I, pos, idx)...]
end

@inline function sample(::LinearResampler, data, I, pos, n::NTuple{Ns,Int}, p::SVector{Ns,T},
                        fillvalue::V) where {Ns,T,V}
  inside(p, n) || return fillvalue
  f = ntuple(k -> floor(p[k]), Val(Ns))
  # in the outer half voxel the neighbor beyond the edge is clamped to the edge voxel
  i0 = ntuple(k -> clamp(unsafe_trunc(Int, f[k]), 1, n[k]), Val(Ns))
  i1 = ntuple(k -> clamp(unsafe_trunc(Int, f[k]) + 1, 1, n[k]), Val(Ns))
  w = ntuple(k -> p[k] - f[k], Val(Ns))
  acc = zero(V)
  for c in 0:(2^Ns-1)
    idx = ntuple(k -> ((c >> (k - 1)) & 1) == 1 ? i1[k] : i0[k], Val(Ns))
    wc = prod(ntuple(k -> ((c >> (k - 1)) & 1) == 1 ? w[k] : one(T) - w[k], Val(Ns)))
    acc += wc * @inbounds(data[with_grid(I, pos, idx)...])
  end
  acc
end

@kernel inbounds = true function resample_kernel!(out, @Const(data), imap, r, pos, n, fillvalue, sd)
  I = @index(Global, Cartesian)
  It = Tuple(I)
  p = imap(grid_index(It, sd))
  out[I] = sample(r, data, It, pos, n, p, fillvalue)
end

function resample!(r::Union{NearestResampler,LinearResampler}, out::AbstractArray, img::AbstractTomoImage,
                   imap::IndexMap; fillvalue=zero(eltype(out)))
  data = parent(img)
  sd = griddims(img)
  pos = ntuple(d -> something(findfirst(==(d), sd), 0), ndims(data))
  n = map(d -> size(data, d), sd)
  backend = KA.get_backend(out)
  resample_kernel!(backend)(out, data, imap, r, Val(pos), n, convert(eltype(out), fillvalue), Val(sd);
                            ndrange=size(out))
  out
end
