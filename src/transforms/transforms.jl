"""
    AbstractSpatialTransform

Supertype of transforms: maps from positions in one space (its [`source`](@ref)) to positions in
another (its [`target`](@ref)). A transform implements `source(t)`, `target(t)` and `t(x)` for a
point `x`, and `inv(t)` where it has an inverse. The functions of this package that take transforms
([`resample`](@ref), [`TomoImages.compose`](@ref), [`transform_geometry`](@ref)) accept
[`AffineTransform`](@ref)s.

Registration results are stored in the direction resampling evaluates: from the space of the fixed
image (the output grid) to the space of the moving image (the image that is resampled). See
[`resample`](@ref).
"""
abstract type AbstractSpatialTransform end

"""
    AffineTransform(A, b; source, target)

The transform `x -> A * x + b` from positions in the space `source` to positions in the space
`target`. The spaces are part of the type, so that [`TomoImages.compose`](@ref) and
[`resample`](@ref) can check them; with `isbits` spaces the transform is `isbits`. `inv(t)` is the
transform from `target` back to `source`.

# Example
```julia
# a rotation by 90 degrees about z and a shift of 120 mm along z, from LPS into the space of a scanner
t = AffineTransform([0 -1 0; 1 0 0; 0 0 1], [0, 0, 120]; source=LPS(), target=MPIScanner())
t([1.0, 0.0, 0.0])                          # [0, 1, 120]
```
"""
struct AffineTransform{W,T<:AbstractFloat,S,D,L} <: AbstractSpatialTransform
  A::SMatrix{W,W,T,L}
  b::SVector{W,T}
  source::S
  target::D
end

function AffineTransform(A::AbstractMatrix, b::AbstractVector; source, target)
  W = length(b)
  size(A) == (W, W) || throw(DimensionMismatch("A must be $W x $W for a shift of length $W, got $(size(A))"))
  T = float(promote_type(eltype(A), eltype(b)))
  AffineTransform{W,T,typeof(source),typeof(target),W * W}(SMatrix{W,W,T}(A), SVector{W,T}(b), source, target)
end

"""
    source(t)

The space a transform maps from.
"""
source(t::AffineTransform) = t.source

"""
    target(t)

The space a transform maps to.
"""
target(t::AffineTransform) = t.target

(t::AffineTransform{W,T})(x) where {W,T} = t.A * SVector{W,T}(Tuple(x)) + t.b

function Base.inv(t::AffineTransform{W,T}) where {W,T}
  Ai = inv(t.A)
  AffineTransform{W,T,typeof(t.target),typeof(t.source),W * W}(Ai, -(Ai * t.b), t.target, t.source)
end

"""
    checkspaces(expected, actual, what)

Throws an `ArgumentError` naming `what` if the spaces `expected` and `actual` differ.
"""
function checkspaces(expected, actual, what)
  expected == actual ||
    throw(ArgumentError("$what: expected space $(spacename(expected)), got $(spacename(actual))"))
  nothing
end

"""
    compose(t2, t1)

    t2 ∘ t1

The transform that applies `t1` and then `t2`: from `source(t1)` to `target(t2)`. Requires
`target(t1) == source(t2)`.
"""
function compose(t2::AffineTransform{W}, t1::AffineTransform{W}) where {W}
  checkspaces(target(t1), source(t2), "compose(t2, t1) needs the source of t2 to be the target of t1")
  A = t2.A * t1.A
  b = t2.A * t1.b + t2.b
  T = eltype(A)
  AffineTransform{W,T,typeof(t1.source),typeof(t2.target),W * W}(A, b, t1.source, t2.target)
end

# Base's composition operator, written without the Unicode character
(::typeof(getfield(Base, Symbol("\u2218"))))(t2::AffineTransform, t1::AffineTransform) = compose(t2, t1)

"""
    transform_between(from, to)

The transform between two spaces whose relation is fixed by convention: [`LPS`](@ref) to
[`RAS`](@ref) and back. Packages add methods for their own spaces, with acquisition metadata as an extra
argument where the relation depends on it: `transform_between(from, to, meta)`.
"""
function transform_between end

transform_between(s::LPS, t::RAS) = AffineTransform(Diagonal(SVector(-1.0, -1.0, 1.0)), zeros(3); source=s, target=t)
transform_between(s::RAS, t::LPS) = AffineTransform(Diagonal(SVector(-1.0, -1.0, 1.0)), zeros(3); source=s, target=t)

"""
    transform_geometry(g, t)

The geometry of the same voxels as `g`, written in the space `target(t)`: the origin is mapped by
`t` and the direction by its linear part, which must be orthogonal (a rotation or a reflection).
The voxels do not move; only the coordinates they are described in change. The result has the
element type of the transform. Throws an `ArgumentError` if the spaces do not fit or the linear
part is not orthogonal (to within `64 * eps(Float32)`).
"""
function transform_geometry(g::GridGeometry{N,W,T}, t::AffineTransform{W}) where {N,W,T}
  checkspaces(space(g), source(t), "transform_geometry needs a transform from the space of the grid")
  A = SMatrix{W,W,Float64}(t.A)
  maximum(abs, A' * A - I) <= 64 * eps(Float32) ||
    throw(ArgumentError("transform_geometry needs a rotation or reflection; the linear part of the transform " *
                        "is $(Matrix(t.A))"))
  GridGeometry(size(g); spacing=Tuple(spacing(g)), origin=t(origin(g)), direction=t.A * direction(g), space=target(t))
end
