"""
    AbstractAxis

Supertype of the non-spatial axes of an image (time, contrast, coil, echo). An axis implements
`length`, `getindex` with a range (the axis of a selection), [`TomoImages.indexof`](@ref) and either
[`TomoImages.coords`](@ref) and [`TomoImages.unit`](@ref) (ordered axes) or
[`TomoImages.labels`](@ref) (unordered axes).
"""
abstract type AbstractAxis end

"""
    NumericAxis(values; unit=nothing)

An ordered axis with a coordinate for every index, in its own unit: time in s, echo time in ms,
frequency in kHz. The coordinates never mix with those of other axes; `unit` is any value that
describes them (a string, or a unit of Unitful.jl).
"""
struct NumericAxis{V<:AbstractVector,U} <: AbstractAxis
  values::V
  unit::U
end
NumericAxis(values::AbstractVector; unit=nothing) = NumericAxis(values, unit)

"""
    CategoricalAxis(labels)

An unordered axis of labels without distances between them: contrasts or tracers, coils.
"""
struct CategoricalAxis{V<:AbstractVector} <: AbstractAxis
  labels::V
end

Base.length(a::NumericAxis) = length(a.values)
Base.length(a::CategoricalAxis) = length(a.labels)
Base.getindex(a::NumericAxis, r::AbstractRange{<:Integer}) = NumericAxis(a.values[r], a.unit)
Base.getindex(a::CategoricalAxis, r::AbstractRange{<:Integer}) = CategoricalAxis(a.labels[r])

"""
    coords(axis)

The coordinates of a [`NumericAxis`](@ref).
"""
coords(a::NumericAxis) = a.values

"""
    unit(axis)

The unit of the coordinates of a [`NumericAxis`](@ref); `nothing` for axes without one.
"""
unit(a::NumericAxis) = a.unit
unit(::CategoricalAxis) = nothing

"""
    labels(axis)

The labels of a [`CategoricalAxis`](@ref).
"""
labels(a::CategoricalAxis) = a.labels

"""
    indexof(axis, value)

The index of `value` on `axis`. For a [`NumericAxis`](@ref) the index of the nearest coordinate,
which must lie within half the distance to its neighbors; for a [`CategoricalAxis`](@ref) the index
of the equal label.
"""
function indexof(a::NumericAxis, v)
  x = a.values
  k = argmin(abs.(x .- v))
  if length(x) == 1
    x[1] == v || throw(ArgumentError("$v is not the coordinate $(x[1]) of the axis"))
    return k
  end
  # inside the axis the nearest coordinate is always within half a gap; beyond the ends it may not be
  gap = k == 1 ? abs(x[2] - x[1]) : abs(x[k] - x[k-1])
  abs(x[k] - v) <= gap / 2 ||
    throw(ArgumentError("$v lies outside the coordinates $(first(x)) to $(last(x)) of the axis"))
  k
end

function indexof(a::CategoricalAxis, label)
  k = findfirst(isequal(label), a.labels)
  k === nothing && throw(ArgumentError("no label $(repr(label)) on the axis, labels: $(a.labels)"))
  k
end
