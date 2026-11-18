# # Define a Space
#
# A space names the coordinate system of a grid or the ends of a transform. TomoImages defines the
# patient spaces [`LPS`](@ref) and [`RAS`](@ref) and a [`DefaultSpace`](@ref); a scanner, a
# session or an atlas gets its own space in the package that knows it. This guide defines spaces,
# gives them names, and adds the transforms that are known from acquisition metadata.
# [Grids, Spaces and Transforms](@ref) explains what spaces are for.

using TomoImages

# ## A space for a scanner
#
# A space is a value of its own type, usually a subtype of [`AbstractSpace`](@ref). A type without
# fields stands for one coordinate system:

struct Scanner <: AbstractSpace end

g = GridGeometry((64, 64, 32); spacing=(0.5, 0.5, 1.0), origin=(-15.75, -15.75, -15.5), space=Scanner())
space(g)

# Keep spaces `isbits` (no strings or arrays in their fields): geometries and transforms then stay
# `isbits` and can be passed into your own GPU kernels. [`resample`](@ref) passes only its index map
# to kernels, so it works with any space.

isbits(g)

# ## Spaces with values
#
# Two sessions on the same scanner are two coordinate systems if the patient moves in between. A
# field distinguishes them; spaces compare by value:

struct Session <: AbstractSpace
  id::Int
end
Session(1) == Session(1), Session(1) == Session(2)

# [`TomoImages.spacename`](@ref) names a space in messages and displays; by default the name of its
# type. A method for your own type gives a better one:

TomoImages.spacename(s::Session) = "session $(s.id)"
TomoImages.spacename(Session(2))

# ## Transforms known from metadata
#
# Some relations between spaces follow from the acquisition, such as the scanner coordinates in the
# patient space from the position of the bed. A method of [`transform_between`](@ref) with the
# metadata as third argument makes it available. Its arguments are types of your package, so it
# extends TomoImages without defining methods for types of other packages:

struct BedPosition
  z::Float64       # the table position along the scanner axis, in mm
end

# For a scanner whose axes are parallel to those of the patient space, the relation is a shift along
# `z` by the bed position (other orientations of the patient add a rotation or reflection):

TomoImages.transform_between(s::Scanner, p::LPS, bed::BedPosition) =
  AffineTransform([1 0 0; 0 1 0; 0 0 1], [0.0, 0.0, bed.z]; source=s, target=p)

t = transform_between(Scanner(), LPS(), BedPosition(120.0))

# [`transform_geometry`](@ref) writes the grid in the patient space; its voxels stay where they are,
# also the far corner:

h = transform_geometry(g, t)
space(h), index_to_world(h, (64, 64, 32)) ≈ t(index_to_world(g, (64, 64, 32)))

# Transforms chain into other spaces: `transform_between(LPS(), RAS()) ∘ t` maps scanner coordinates
# to RAS. To resample an image of the scanner onto a grid in LPS, [`resample`](@ref) needs the
# transform from the grid's space to the image's space, here `inv(t)` (from LPS to the scanner).
#
# [`resample`](@ref), [`transform_geometry`](@ref) and [`TomoImages.compose`](@ref) check spaces.
# Composing transforms whose spaces do not chain is an error:

try #hide
TomoImages.compose(t, t)
catch e #hide
  print(sprint(showerror, e)) #hide
end #hide
