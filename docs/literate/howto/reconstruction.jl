# # Wrap a Reconstruction in a SpatialImage
#
# A reconstruction ends with an array. This guide, for developers of reconstruction packages, turns
# it into a [`SpatialImage`](@ref) that carries where its voxels lie, what its dimensions are and
# how it was made: an MRI reconstruction of an oblique slab with echoes and coils, and an MPI
# reconstruction with contrasts and a time series.

using TomoImages

# ## What a scanner reports
#
# Scanner metadata usually describes the reconstructed volume by its field of view, its matrix
# size, the position of its center and the directions of its axes. A [`GridGeometry`](@ref) needs
# the spacing, the center of the *first* voxel and the directions:
#
# * spacing = field of view ./ matrix size (the field of view spans from the outer edge of the first
#   voxel to the outer edge of the last);
# * origin = center - direction * (spacing .* (matrix size .- 1) ./ 2).
#
# ## MRI: an oblique slab
#
# An MRI reconstruction of 64 × 64 × 8 voxels over a field of view of 220 × 220 × 40 mm, three
# echoes and four receive coils. The slab is rotated by 20 degrees about the `z` axis of the patient
# space [`LPS`](@ref); its read, phase and slice directions are the columns of `D`:

matrix = (64, 64, 8)
fov = (220.0, 220.0, 40.0)
center = [10.0, -5.0, 30.0]
dir_read, dir_phase, dir_slice = [cosd(20), sind(20), 0.0], [-sind(20), cosd(20), 0.0], [0.0, 0.0, 1.0]
D = hcat(dir_read, dir_phase, dir_slice)
voxel = fov ./ matrix
g = GridGeometry(matrix; spacing=voxel, origin=center - D * collect(voxel .* (matrix .- 1) ./ 2), direction=D, space=LPS())

# The center of the grid is where the scanner put it:

fov_center(g) ≈ center

# The reconstruction returns one complex value per voxel, echo and coil. The dimensions get the
# names of the MRI axes; the echo dimension holds the echo times, the coil dimension the coil
# numbers as labels:

reco = zeros(ComplexF32, 64, 64, 8, 3, 4)
img = SpatialImage(reco, g;
                   dims=(read=GridAxis(1), phase=GridAxis(2), slice=GridAxis(3),
                         echo=NumericAxis([5.0, 10.0, 15.0]; unit="ms"), coil=CategoricalAxis(1:4)),
                   meta=(quantity="MR signal", unit="a.u.", method="direct reconstruction"))
summary(img)

# A user of the image selects by name or by value and gets correctly placed images:

summary(selectvalue(img; echo=10.0, coil=2)), space(img[slice=4])

# ## MPI: contrasts and a time series
#
# An MPI reconstruction is placed in the coordinates of its scanner. The package that reads the
# scanner data defines a space for it, a type of its own (see [Define a Space](@ref)):

struct MPIScanner <: AbstractSpace end

# 20 × 20 × 10 voxels over 40 × 40 × 20 mm, centered 5 mm above the center of the scanner along
# `z`; two tracers and 50 frames of 21.5 ms. The grid dimensions sit between the contrast and the
# time dimension, in the order the reconstruction produces them:

matrix = (20, 20, 10)
voxel = (40.0, 40.0, 20.0) ./ matrix
g = GridGeometry(matrix; spacing=voxel, origin=[0.0, 0.0, 5.0] .- voxel .* (matrix .- 1) ./ 2, space=MPIScanner())
c = rand(Float32, 2, 20, 20, 10, 50)
mpi = SpatialImage(c, g;
                   dims=(contrast=CategoricalAxis(["tracer A", "tracer B"]),
                         x=GridAxis(1), y=GridAxis(2), z=GridAxis(3),
                         time=NumericAxis(0.0215 .* (0:49); unit="s")),
                   meta=(quantity="concentration", unit="mmol/L"))
summary(mpi)

# The same functions apply as for the MRI image:

summary(selectvalue(mpi; contrast="tracer B", time=0.5)), TomoImages.spatialdims(mpi)

# ## Checklist
#
# * The origin is the center of the first voxel. Convert from a center of the field of view as
#   above.
# * Spacings are positive; a flipped axis is a direction column with the opposite sign.
# * The columns of the direction are orthonormal. Values read from files in `Float32` are
#   accepted.
# * The space says which coordinate system the positions are in: [`LPS`](@ref) or
#   [`RAS`](@ref) where the scanner reports patient coordinates, a space of the package for scanner
#   coordinates. Where the metadata relates the scanner to the patient, a method of
#   [`transform_between`](@ref) makes that available.
# * Every dimension has a name; non-grid dimensions have coordinates with a unit
#   ([`NumericAxis`](@ref)) or labels ([`CategoricalAxis`](@ref)).
# * The metadata says what the values are and how they were computed.
