# The axes of an image: one entry per array dimension. Grid axes refer to an index axis of the
# geometry; every other dimension is described by a vector with one value per index.

"""
    GridAxis(k)

Marks an array dimension as index axis `k` of the image's geometry. Grid axes carry no coordinates
of their own: the geometry places all of them together, since its direction may mix them. The name
of the dimension (the key in `dims`) only labels the index axis; it does not say that the axis runs
along an axis of the space.
"""
struct GridAxis
  index::Int
end

"""
    NumericAxis(values; unit=nothing)

An ordered axis with a coordinate for every index, in its own unit: time in s, echo time in ms,
frequency in kHz. It is an `AbstractVector` of its coordinates, which should be sorted (selection by
value assumes it; it is not checked). `unit` is any value that describes
them (a string, or a unit of Unitful.jl).

Any `AbstractVector` can describe a non-grid dimension; `NumericAxis` adds the unit.
"""
struct NumericAxis{T,V<:AbstractVector{T},U} <: AbstractVector{T}
  values::V
  unit::U
end
NumericAxis(values::AbstractVector; unit=nothing) = NumericAxis(values, unit)

"""
    CategoricalAxis(labels)

An unordered axis of labels without distances between them: contrasts or tracers, coils. It is an
`AbstractVector` of its labels. Labels are found by equality, also when they are numbers (coil
`3`).
"""
struct CategoricalAxis{T,V<:AbstractVector{T}} <: AbstractVector{T}
  labels::V
end

Base.size(a::NumericAxis) = size(a.values)
Base.size(a::CategoricalAxis) = size(a.labels)
Base.@propagate_inbounds Base.getindex(a::NumericAxis, i::Int) = a.values[i]
Base.@propagate_inbounds Base.getindex(a::CategoricalAxis, i::Int) = a.labels[i]
Base.getindex(a::NumericAxis, r::AbstractRange{<:Integer}) = NumericAxis(a.values[r], a.unit)
Base.getindex(a::CategoricalAxis, r::AbstractRange{<:Integer}) = CategoricalAxis(a.labels[r])

"""
    unit(axis)

The unit of the coordinates of a [`NumericAxis`](@ref); `nothing` for other axes.
"""
unit(a::NumericAxis) = a.unit
unit(::AbstractVector) = nothing

"""
    indexof(axis, value)

The index of `value` on the axis. On a numeric axis (a vector of numbers that is not a
[`CategoricalAxis`](@ref)) the index of the nearest coordinate, which must lie within half the
distance to its neighbor at the ends (on an axis of one coordinate, the value must be equal); on
other axes the index of the first equal element. Throws an `ArgumentError` if there is none.
"""
function indexof(a::AbstractVector{<:Number}, v::Number)
  k = argmin(abs.(a .- v))
  if length(a) == 1
    a[k] == v || throw(ArgumentError("$v is not the coordinate $(a[k]) of the axis"))
    return k
  end
  # inside the axis the nearest coordinate is always within half a gap; beyond the ends it may not be
  gap = k == firstindex(a) ? abs(a[k+1] - a[k]) : abs(a[k] - a[k-1])
  abs(a[k] - v) <= gap / 2 ||
    throw(ArgumentError("$v lies outside the coordinates $(first(a)) to $(last(a)) of the axis"))
  k
end
indexof(a::CategoricalAxis, v) = findlabel(a, v)
# numbers on a categorical axis are labels
indexof(a::CategoricalAxis{<:Number}, v::Number) = findlabel(a, v)
indexof(a::AbstractVector, v) = findlabel(a, v)

function findlabel(a, v)
  k = findfirst(isequal(v), a)
  k === nothing && throw(ArgumentError("no element $(repr(v)) on the axis $(collect(a))"))
  k
end

# default names of the dimensions: grid axes x1, x2, ...; other dimensions dim_d
default_dims(::Val{Ns}, N) where {Ns} =
  NamedTuple{ntuple(d -> d <= Ns ? Symbol(:x, d) : Symbol(:dim_, d), N)}(
    ntuple(d -> d <= Ns ? GridAxis(d) : nothing, N))
