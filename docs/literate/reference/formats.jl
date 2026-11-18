# # Conventions of Formats and Packages
#
# File formats and other packages describe the placement of a grid with other numbers. This page
# lists the conventions of the common ones and translates example values of each into a
# [`GridGeometry`](@ref), checked against the format's own formula for the position of a voxel
# (see the references at the end).

using TomoImages, LinearAlgebra

# ## Overview
#
# | | Origin | Index of the first voxel | Orientation | Coordinate system |
# |:---|:---|:---|:---|:---|
# | TomoImages | center of the first voxel | 1 | direction matrix, columns per index axis | space type |
# | DICOM | Image Position (Patient): center of the first voxel | 0 | Image Orientation (Patient): direction cosines of a row and a column | LPS |
# | NIfTI | sform: position of voxel `(0, 0, 0)` | 0 | the sform matrix (spacing included) | RAS |
# | ImageGeoms.jl | offsets, in voxels, of the grid center | 1 | none; spacings may be negative | none |
#
# ## DICOM
#
# A DICOM slice gives the position of the center of its first pixel (Image Position (Patient)), the
# direction cosines of its rows and columns (Image Orientation (Patient)) and the pixel spacing as
# the distance between adjacent rows, then between adjacent columns (Pixel Spacing), in mm in the
# LPS patient space. The position of the pixel in column `i` and row `j`, both counted from 0, is
#
# ```math
# P = S + X \Delta_i \, i + Y \Delta_j \, j,
# ```
#
# with `S` the Image Position, `X` the direction cosines of a row (the direction in which the column
# index grows), `Y` those of a column, `Δi` the column spacing (the second value of Pixel Spacing)
# and `Δj` the row spacing (the first value). Example values of an oblique slice:

image_position = [-120.0, -95.5, 40.0]
row_cosines = [1.0, 0.0, 0.0]
column_cosines = [0.0, cosd(15), -sind(15)]
pixel_spacing = [0.8, 0.7];             # row spacing, column spacing

# For an array indexed (column, row), index axis 1 runs along a row and index axis 2 along a column.
# The slice is a 2D grid in the 3D patient space:

g = GridGeometry((256, 192); spacing=(pixel_spacing[2], pixel_spacing[1]), origin=image_position,
                 direction=hcat(row_cosines, column_cosines), space=LPS())
dicom(i, j) = image_position + row_cosines * pixel_spacing[2] * i + column_cosines * pixel_spacing[1] * j
all(index_to_world(g, (i + 1, j + 1)) ≈ dicom(i, j) for i in (0, 17, 255), j in (0, 101, 191))

# Taking the two values of Pixel Spacing in the other order changes nothing for square pixels and
# stretches the image otherwise.
#
# ## NIfTI
#
# A NIfTI file maps voxel indices `(i, j, k)`, counted from 0, to positions `(x, y, z)` in mm in
# the RAS patient space, through the 4 × 4 matrix of its sform. The first three columns hold the
# steps along the index axes, spacing included, and the fourth the position of voxel `(0, 0, 0)`:

R = [1 0 0; 0 cosd(6) -sind(6); 0 sind(6) cosd(6)]     # a rotation by 6 degrees about x
A = [R * Diagonal([-0.9, 0.9, 0.9]) [90.0, -126.0, -72.0]
     0 0 0 1];

# The spacing is the length of each column, the direction the columns divided by it, the origin the
# last column:

spacing_nifti = [norm(A[1:3, d]) for d in 1:3]
g = GridGeometry((200, 240, 180); spacing=spacing_nifti, origin=A[1:3, 4],
                 direction=A[1:3, 1:3] ./ spacing_nifti', space=RAS())
all(index_to_world(g, (i + 1, j + 1, k + 1)) ≈ (A * [i, j, k, 1])[1:3] for i in (0, 99), j in (0, 239), k in (0, 50))

# The first index axis runs toward `-x` (a negative first entry), which the direction matrix
# expresses with a flipped column; the spacing stays positive. In LPS the same grid has its origin
# at:

TomoImages.origin(transform_geometry(g, transform_between(RAS(), LPS())))

# An sform may also contain a shear, which a [`GridGeometry`](@ref) does not represent: its
# direction must have orthonormal columns.
#
# ## ImageGeoms.jl
#
# [ImageGeoms.jl](https://github.com/JuliaImageRecon/ImageGeoms.jl), with the conventions of the
# Michigan Image Reconstruction Toolbox, places voxel `i` (counted from 1) of a dimension with `n`
# voxels, spacing `Δ` and offset `o` at
#
# ```math
# x = \left((i - 1) - \left(\tfrac{n - 1}{2} + o\right)\right) \Delta,
# ```
#
# so the grid is centered at zero for `o = 0`, and the offset moves it by whole or fractional
# voxels. A negative `Δ` reverses the axis. For a 3D grid with a negative spacing along `y`:

sizes, deltas, offsets = (128, 96, 40), (0.5, -0.5, 1.0), (0.0, 0.5, -2.0)
imagegeoms(i, d) = ((i - 1) - ((sizes[d] - 1) / 2 + offsets[d])) * deltas[d];

# The origin is the position of voxel 1 along every axis; the spacing is the magnitude of `Δ`, and
# its sign goes into the direction:

g = GridGeometry(sizes; spacing=abs.(deltas), origin=[imagegeoms(1, d) for d in 1:3],
                 direction=Diagonal(collect(sign.(deltas))))
all(index_to_world(g, (i, j, k)) ≈ [imagegeoms(i, 1), imagegeoms(j, 2), imagegeoms(k, 3)]
    for i in (1, 128), j in (1, 50), k in (1, 40))

# ImageGeoms.jl has no coordinate system; the geometry gets the [`DefaultSpace`](@ref) unless the
# scanner's space is known.
#
#
# ## References
#
# - DICOM Standard, PS3.3 Information Object Definitions, Section C.7.6.2 Image Plane Module, in
#   particular C.7.6.2.1.1 Image Position and Image Orientation, National Electrical Manufacturers
#   Association:
#   [dicom.nema.org](https://dicom.nema.org/medical/dicom/current/output/chtml/part03/sect_C.7.6.2.html).
# - NIfTI-1 Data Format, header file `nifti1.h` ("Official definition of the nifti1 header"),
#   Neuroimaging Informatics Technology Initiative:
#   [NIFTI-Imaging/nifti_clib](https://github.com/NIFTI-Imaging/nifti_clib/blob/master/niftilib/nifti1.h).
# - ImageGeoms.jl, `src/core.jl` (the function `axis`):
#   [JuliaImageRecon/ImageGeoms.jl](https://github.com/JuliaImageRecon/ImageGeoms.jl).
