"""
    AbstractSpace

Supertype of coordinate spaces: the coordinate system in which the origin and direction of a grid,
and the points a transform maps, are given. A space holds no numbers about other spaces; how two
spaces relate is the job of transforms ([`AffineTransform`](@ref)).

A space is any value with
* `==`: whether two spaces are the same coordinate system (by default, structs compare their
  fields, so a parametric space such as `struct Session <: AbstractSpace; id::Int; end` compares by
  value);
* [`spacename`](@ref): a name for display and files.

Spaces should be `isbits`, so that geometries and transforms holding them can be passed into
kernels.

# Example
```julia
struct MPIScanner <: AbstractSpace end      # defined by the package that owns the space
g = GridGeometry((40, 40, 20); spacing=(1.0, 1.0, 1.0), space=MPIScanner())
```
"""
abstract type AbstractSpace end

"""
    DefaultSpace()

The space of geometries that do not name one: for code that works in a single coordinate system,
where no image ever meets an image from another scanner or session.
"""
struct DefaultSpace <: AbstractSpace end

"""
    LPS()

The patient space of DICOM: `+x` toward the patient's left, `+y` toward posterior (the
back), `+z` toward superior (the head). See [`RAS`](@ref) for the other common convention and
[`transform_between`](@ref) for the transform between both.
"""
struct LPS <: AbstractSpace end

"""
    RAS()

The patient space of NIfTI and most neuroimaging tools: `+x` toward the patient's right, `+y` toward
anterior, `+z` toward superior. A point with coordinates `(x, y, z)` in [`LPS`](@ref) has the
coordinates `(-x, -y, z)` in `RAS`.
"""
struct RAS <: AbstractSpace end

"""
    spacename(space) -> String

The name of `space` for display and files. Defaults to the name of its type.
"""
spacename(s) = string(nameof(typeof(s)))
