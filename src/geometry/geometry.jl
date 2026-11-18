"""
    AbstractGeometry

Supertype of geometries: maps between voxel indices and positions in a coordinate space. A geometry
implements `size`, [`index_to_world`](@ref), [`world_to_index`](@ref) and [`space`](@ref); functions of
this package use only these, never the fields.
"""
abstract type AbstractGeometry end

"""
    GridGeometry(size; spacing, origin, direction, space=DefaultSpace())

A regular grid of `size` voxels (`N` index dimensions) in a coordinate space with `W >= N` world
dimensions. Voxel `i` (1-based) has its center at

    x = origin + direction * (spacing .* (i .- 1))

* `spacing`: the distance between voxel centers along each index axis, positive (default 1).
* `origin`: the center of the first voxel, in the space (default 0); its length is `W`.
* `direction`: a `W × N` matrix whose columns are the unit directions of the index axes in the
  space (default: the first `N` axes of the space). The columns must be orthonormal, so a grid is
  rotated or reflected but never sheared; a flipped axis is a column with the opposite sign, never
  a negative spacing.
* `space`: the coordinate space the origin and direction are given in (an [`AbstractSpace`](@ref)).

`W > N` describes a grid of lower dimension in a larger space, such as a slice of a 3D volume (see
indexing the geometry, [`TomoImages.subgrid`](@ref)). Lengths have no unit attached; all geometries of one application use the same
(mm by convention).

The orthonormality check accepts direction matrices stored in `Float32`: their columns are
orthonormal only up to rounding, a few `eps(Float32)`, which the check tolerates up to
`64 * eps(Float32)`.

# Example
```julia
g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=LPS())
index_to_world(g, (1, 1, 1))                # the origin
fov_center(g)                               # (0, 0, 0)
```
"""
struct GridGeometry{N,W,T<:AbstractFloat,S,L} <: AbstractGeometry
  size::NTuple{N,Int}
  spacing::SVector{N,T}
  origin::SVector{W,T}
  direction::SMatrix{W,N,T,L}
  space::S

  function GridGeometry{N,W,T,S,L}(size, spacing, origin, direction, space) where {N,W,T,S,L}
    W >= N || throw(ArgumentError("a grid of $N index dimensions needs a space of at least $N dimensions, got $W"))
    all(>(0), size) || throw(ArgumentError("size must be positive, got $size"))
    all(s -> s > 0 && isfinite(s), spacing) ||
      throw(ArgumentError("spacing must be positive and finite, got $(Tuple(spacing)); a flipped axis is a " *
                          "direction column with the opposite sign"))
    all(isfinite, origin) || throw(ArgumentError("origin must be finite, got $(Tuple(origin))"))
    D = SMatrix{W,N,Float64}(direction)
    maximum(abs, D' * D - I) <= 64 * eps(Float32) ||
      throw(ArgumentError("the columns of direction must be orthonormal (rotations and reflections only), got $direction"))
    new{N,W,T,S,L}(size, spacing, origin, direction, space)
  end
end

function GridGeometry(size::NTuple{N,Integer}; spacing=ntuple(_ -> 1.0, N), origin=nothing, direction=nothing,
                      space=DefaultSpace()) where {N}
  W = origin !== nothing ? length(origin) : direction !== nothing ? Base.size(direction, 1) : N
  T = float(promote_type(eltype(spacing), origin === nothing ? Float64 : eltype(origin),
                         direction === nothing ? Float64 : eltype(direction)))
  T = T <: AbstractFloat ? T : Float64
  length(spacing) == N || throw(ArgumentError("spacing needs $N entries, got $(length(spacing))"))
  o = origin === nothing ? zeros(SVector{W,T}) : SVector{W,T}(Tuple(origin))
  D = direction === nothing ? SMatrix{W,N,T}(ntuple(k -> T((k - 1) % W == div(k - 1, W)), W * N)) :
      SMatrix{W,N,T}(direction)
  GridGeometry{N,W,T,typeof(space),W * N}(Int.(size), SVector{N,T}(Tuple(spacing)), o, D, space)
end

Base.size(g::GridGeometry) = g.size
Base.size(g::GridGeometry, d::Integer) = g.size[d]
Base.ndims(::GridGeometry{N}) where {N} = N
Base.eltype(::GridGeometry{N,W,T}) where {N,W,T} = T

"""
    worlddims(g)

The number of dimensions of the space of `g`: `ndims(g)` for a full grid, more for a slice.
"""
worlddims(::GridGeometry{N,W}) where {N,W} = W

"""
    spacing(g)

The distance between voxel centers along each index axis of `g`.
"""
spacing(g::GridGeometry) = g.spacing

"""
    origin(g)

The center of the first voxel of `g`, in its space.
"""
origin(g::GridGeometry) = g.origin

"""
    direction(g)

The `W × N` matrix whose columns are the unit directions of the index axes of `g` in its space.
"""
direction(g::GridGeometry) = g.direction

"""
    space(g)

The coordinate space of a geometry, or of the geometry of an image.
"""
space(g::GridGeometry) = g.space

"""
    index_to_world(g, i)

The position of the (possibly fractional) voxel index `i` of `g` in its space, as an `SVector`.
Integer indices give voxel centers.
"""
@inline index_to_world(g::GridGeometry{N,W,T}, i) where {N,W,T} =
  g.origin + g.direction * (g.spacing .* (SVector{N,T}(Tuple(i)) .- one(T)))

"""
    world_to_index(g, x)

The continuous voxel index of the position `x` in the space of `g`: the inverse of
[`index_to_world`](@ref). For a grid of lower dimension than its space (a slice), `x` is first
projected onto the plane of the grid.
"""
@inline world_to_index(g::GridGeometry{N,W,T}, x) where {N,W,T} =
  (g.direction' * (SVector{W,T}(Tuple(x)) - g.origin)) ./ g.spacing .+ one(T)

"""
    affine(g)

The `(W + 1) × (N + 1)` matrix of [`index_to_world`](@ref) in homogeneous coordinates: `affine(g) *
[i; 1] == [x; 1]` for the 1-based index `i` and the position `x`.
"""
function affine(g::GridGeometry{N,W,T}) where {N,W,T}
  M = g.direction * Diagonal(g.spacing)
  c = g.origin - M * ones(SVector{N,T})
  top = hcat(M, c)
  bottom = hcat(zeros(SMatrix{1,N,T}), SMatrix{1,1,T}(1))
  vcat(top, bottom)
end

"""
    fov(g)

The extent of the grid along each index axis, from the outer edge of the first voxel to the outer
edge of the last: `spacing .* size`.
"""
fov(g::GridGeometry) = g.spacing .* SVector(g.size)

"""
    fov_center(g)

The position of the center of the grid in its space.
"""
fov_center(g::GridGeometry{N,W,T}) where {N,W,T} = index_to_world(g, (SVector{N,T}(g.size) .+ one(T)) ./ 2)

"""
    axiscoords(g, d)

The coordinates of the voxel centers along index axis `d`, as a range, for an axis that runs along
axis `d` of the space (column `d` of the direction is the unit vector `e_d`). Packages that take
coordinate vectors per axis, such as ImagePhantoms.jl, use it for axis-aligned grids.
"""
function axiscoords(g::GridGeometry{N,W,T}, d::Integer) where {N,W,T}
  e = SVector{W,T}(ntuple(k -> T(k == d), W))
  g.direction[:, d] == e ||
    throw(ArgumentError("index axis $d of the grid does not run along axis $d of its space; use index_to_world"))
  range(g.origin[d]; step=g.spacing[d], length=g.size[d])
end

"""
    isapprox(a::GridGeometry, b::GridGeometry; atol)

Whether `a` and `b` describe the same voxels: the same size and space, and voxel centers that agree
to within `atol` (a length in the unit of the geometries; default a thousandth of the smallest voxel
spacing of `a`). The centers are compared at the corners of the grid, which determine an affine
grid completely.
"""
function Base.isapprox(a::GridGeometry{N,W}, b::GridGeometry; atol=1e-3 * minimum(spacing(a))) where {N,W}
  (size(a) == size(b) && worlddims(b) == W && space(a) == space(b)) || return false
  for c in CartesianIndices(ntuple(_ -> 2, N))
    i = ntuple(d -> c[d] == 1 ? 1 : size(a, d), N)
    sqrt(sum(abs2, index_to_world(a, i) - index_to_world(b, i))) <= atol || return false
  end
  true
end

"""
    in(x, g::GridGeometry)

Whether the position `x` lies inside the voxels of `g`: within half a voxel beyond the outermost
voxel centers along every index axis. For a grid of lower dimension than its space (a slice), `x` is
projected onto the plane of the grid.
"""
function Base.in(x, g::GridGeometry{N}) where {N}
  i = world_to_index(g, x)
  all(k -> 0.5 <= i[k] <= size(g, k) + 0.5, 1:N)
end

Base.CartesianIndices(g::GridGeometry) = CartesianIndices(size(g))

Base.getindex(g::GridGeometry{N}, inds::Vararg{Union{Integer,AbstractRange{<:Integer},Colon},N}) where {N} =
  subgrid(g, ntuple(d -> inds[d] isa Colon ? (1:size(g, d)) : inds[d], N))

"""
    g[inds...]
    subgrid(g, inds)

The geometry of the voxels `inds` of `g`: one entry per index dimension, an integer (the dimension
is dropped), a range (kept, with its step) or a colon. Every selected voxel keeps its position in the space:
the origin moves to the first selected voxel, the spacing is scaled by the step, and a range with a
negative step flips the direction column. Dropping a dimension keeps the space, so the result may be
a slice in a space of more dimensions.
"""
function subgrid(g::GridGeometry{N,W,T}, inds::Tuple) where {N,W,T}
  length(inds) == N || throw(ArgumentError("need $N indices, got $(length(inds))"))
  for (d, r) in enumerate(inds)
    if r isa Integer
      1 <= r <= size(g, d) || throw(BoundsError(g, inds))
    elseif r isa AbstractRange{<:Integer}
      isempty(r) && throw(ArgumentError("an empty range selects no voxels along dimension $d"))
      (1 <= minimum(r) && maximum(r) <= size(g, d)) || throw(BoundsError(g, inds))
    else
      throw(ArgumentError("indices of a grid are integers or integer ranges, got $(typeof(r))"))
    end
  end
  first_index = ntuple(d -> inds[d] isa Integer ? Int(inds[d]) : Int(first(inds[d])), N)
  kept = [d for d in 1:N if inds[d] isa AbstractRange]
  n = length(kept)
  sz = ntuple(k -> length(inds[kept[k]]), n)
  sp = ntuple(k -> g.spacing[kept[k]] * abs(step(inds[kept[k]])), n)
  D = SMatrix{W,n,T}(reduce(hcat, [sign(step(inds[d])) * g.direction[:, d] for d in kept]; init=zeros(T, W, 0)))
  GridGeometry(sz; spacing=sp, origin=index_to_world(g, first_index), direction=D, space=g.space)
end

function Base.show(io::IO, ::MIME"text/plain", g::GridGeometry{N,W,T}) where {N,W,T}
  aligned = g.direction == SMatrix{W,N,T}(ntuple(k -> T((k - 1) % W == div(k - 1, W)), W * N))
  println(io, "GridGeometry ", join(g.size, "x"), " in ", spacename(g.space), " (", W, "D)")
  println(io, "  spacing:   ", Tuple(g.spacing))
  println(io, "  origin:    ", Tuple(g.origin))
  print(io, "  direction: ", aligned ? "axis-aligned" : string(Matrix(g.direction)))
end
