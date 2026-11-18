"""
    IndexMap{N,T}

The affine map from a voxel index `j` of the output grid to the continuous voxel index `i = M * j +
c` of the input grid at the same position (after the transform, if any). `isbits`, so that kernels
take it as an argument; `m(j)` evaluates it.
"""
struct IndexMap{N,T,L}
  M::SMatrix{N,N,T,L}
  c::SVector{N,T}
end

@inline (m::IndexMap{N,T})(j) where {N,T} = m.M * SVector{N,T}(Tuple(j)) + m.c

"""
    index_map(out, in[, t]; eltype=Float64)

The [`IndexMap`](@ref) from the voxels of the geometry `out` to the continuous voxel indices of the
geometry `in`: the position of output voxel `j` in the space of `out`, mapped by `t` into the space
of `in` (without `t`, positions are taken as they are), as an index of `in`. Both geometries need
the same numbers of index and world dimensions. The map is computed in `Float64` and stored in
`eltype`. It does not check spaces; [`resample`](@ref) does before it calls `index_map`.
"""
function index_map(out::GridGeometry{N,W}, in::GridGeometry{N,W}, t=nothing; eltype::Type{T}=Float64) where {N,W,T}
  A, b = transform_parts(t, W)
  So = Diagonal(SVector{N,Float64}(spacing(out)))
  Si = Diagonal(1 ./ SVector{N,Float64}(spacing(in)))
  Do = SMatrix{W,N,Float64}(direction(out))
  Di = SMatrix{W,N,Float64}(direction(in))
  oo = SVector{W,Float64}(origin(out))
  oi = SVector{W,Float64}(origin(in))
  # i = Si Di' (A (oo + Do So (j - 1)) + b - oi) + 1
  M = Si * Di' * A * Do * So
  c = Si * (Di' * (A * (oo - Do * So * ones(SVector{N,Float64})) + b - oi)) .+ 1
  IndexMap{N,T,N * N}(SMatrix{N,N,T}(M), SVector{N,T}(c))
end

transform_parts(::Nothing, W) = SMatrix{W,W,Float64}(I), zeros(SVector{W,Float64})
transform_parts(t::AffineTransform{W}, ::Integer) where {W} = SMatrix{W,W,Float64}(t.A), SVector{W,Float64}(t.b)
