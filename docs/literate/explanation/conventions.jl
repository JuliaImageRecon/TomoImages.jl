# # Grids, Spaces and Transforms
#
# This page explains what the numbers of a grid geometry mean, what a space adds to them, and how
# transforms connect spaces and are used by resampling.

using TomoImages, CairoMakie
include(joinpath(pkgdir(TomoImages), "docs", "utils", "figures.jl")) #hide
nothing #hide

# ## Grids
#
# A [`GridGeometry`](@ref) of `N` index dimensions places voxel `i` (1-based) at
#
# ```math
# x = o + D \, \operatorname{diag}(s) \, (i - 1)
# ```
#
# * the **origin** ``o`` is the center of the first voxel, not a corner and not the center of the
#   grid;
# * the **spacing** ``s`` is the distance between neighboring voxel centers along each index axis,
#   always positive;
# * the **direction** ``D`` holds the unit directions of the index axes as columns. Its columns are
#   orthonormal: a grid can be rotated or reflected, never sheared or scaled.
#
# This is the convention of ITK. A grid whose index axis 1 runs along `-x`
# has a direction column `(-1, 0)`, never a negative spacing.
#
# The figure shows a grid of 6 × 4 pixels of 1 mm × 1.5 mm, rotated by 20 degrees: the gray lines are
# the pixel edges, the dots the pixel centers.

c, s = cos(deg2rad(20)), sin(deg2rad(20))
g = GridGeometry((6, 4); spacing=(1.0, 1.5), origin=(1.0, 0.5), direction=[c -s; s c])
fig = grid_figure(g) #hide
save("conventions_grid.png", fig) #hide
nothing #hide
# ![A rotated grid with its origin and direction](conventions_grid.png)

# Derived quantities are functions, never stored: [`affine`](@ref) (the matrix in homogeneous
# coordinates), [`fov`](@ref) (the extent from edge to edge), [`fov_center`](@ref).
#
# ## Selections and slices
#
# Selecting voxels keeps their positions ([`subgrid`](@ref)): a range moves the origin to its first
# voxel and multiplies the spacing by its step, a reversed range flips the direction column. An
# integer drops the dimension, and the result stays in the space of the full grid: a slice of a
# rotated volume is a 2D grid in a 3D space, whose direction is a 3 × 2 matrix. Its pixels keep
# their 3D positions:

volume = GridGeometry((8, 8, 4); spacing=(1.0, 1.0, 2.0), direction=[c -s 0; s c 0; 0 0 1])
slice = subgrid(volume, (1:8, 1:8, 3))
size(direction(slice)), index_to_world(slice, (2, 5)) == index_to_world(volume, (2, 5, 3))

# ## Spaces
#
# The positions of a grid are numbers, and numbers only mean something in a coordinate system: the
# coordinates of a CT scanner, those of an MPI device, or a patient space such as [`LPS`](@ref)
# (`+x` toward the patient's left, `+y` posterior, `+z` superior) or [`RAS`](@ref) (`+x` right,
# `+y` anterior, `+z` superior). The space of a grid names this coordinate system.
#
# Two grids with the same numbers in different spaces describe different voxels, so
# [`samegrid`](@ref) and [`resample`](@ref) treat them as different. Conversely, the same voxels
# have different numbers in different spaces. In LPS and RAS, the first two coordinates change
# sign:

lps = GridGeometry((4, 4, 2); origin=(10.0, 20.0, 30.0), space=LPS())
ras = transform_geometry(lps, transform_between(LPS(), RAS()))
origin(ras), index_to_world(ras, (4, 1, 1)) == transform_between(LPS(), RAS())(index_to_world(lps, (4, 1, 1)))

# A space holds no numbers itself. Spaces are types, so packages define their own
# ([Define a Space](@ref)), and they are part of the types of geometries and transforms.
#
# ## Transforms
#
# A transform maps positions from its [`source`](@ref) space to its [`target`](@ref) space. The
# relations come from three places:
#
# * conventions that hold without data, such as LPS and RAS ([`transform_between`](@ref));
# * acquisition metadata, such as the position of the bed (methods of `transform_between` in the
#   packages that read the metadata);
# * registration or calibration between devices, computed by other packages.
#
# [`compose`](@ref) chains transforms and checks that the target of the first is the source of the
# second.
#
# ## Resampling pulls
#
# [`resample`](@ref) computes, for every voxel of the output grid, where to read in the input:
# output index, position in the output space, the transform into the input space, input index. Each
# output voxel gets exactly one value. The opposite, sending every input voxel to the output,
# leaves output voxels without a value and gives others several, as soon as the grids differ in
# spacing or orientation.
#
# Resampling therefore needs the transform from the output space to the input space. In
# registration terms the output is the grid of the fixed image and the input is the moving image,
# and registration results are stored in this direction, from fixed to moving. For an affine
# transform this is a convention, since its inverse is cheap. For a deformation field it decides
# whether resampling needs an inverse that has no closed form.
#
# The resamplers read the input over the extent of its voxels: up to half a voxel beyond the
# outermost voxel centers they take the value at the edge, and farther out they write a fill value.
# Nearest-neighbor resampling keeps the element type and is the default for integer images such as
# label maps, whose values must not be averaged. Multilinear interpolation is the default
# otherwise.
