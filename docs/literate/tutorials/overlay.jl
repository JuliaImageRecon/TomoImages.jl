# # Overlay Images from Two Devices
#
# A CT image shows the anatomy, an MPI image the distribution of a tracer. Each device reconstructs
# in its own coordinates, so the two images cannot be overlaid voxel by voxel, even where their
# numbers look alike. This tutorial gives each image its space, relates the spaces by a transform
# and resamples the MPI image onto the CT grid. The example is 2D; the same calls work in 3D.

using TomoImages, CairoMakie
include(joinpath(pkgdir(TomoImages), "docs", "utils", "figures.jl")) #hide
nothing #hide

# ## Two spaces
#
# Each device gets a space: a type, so that images and transforms carry it and functions can check
# it (see [Define a Space](@ref)).

struct CTScanner <: AbstractSpace end
struct MPIDevice <: AbstractSpace end

# ## The CT image
#
# 96 × 96 pixels of 0.5 mm in the CT scanner's coordinates. The anatomy: a body (attenuation 0.2)
# with two bones (1.0).

ct_grid = GridGeometry((96, 96); spacing=(0.5, 0.5), origin=(-23.75, -23.75), space=CTScanner())
body(x) = (x[1] / 20)^2 + (x[2] / 14)^2 < 1 ? 0.2 : 0.0
bones(x) = hypot(x[1] + 8, x[2] - 2) < 3 || hypot(x[1] - 9, x[2] + 1) < 2.5 ? 0.8 : 0.0
ct = SpatialImage([body(x) + bones(x) for x in (index_to_world(ct_grid, (i, j)) for i in 1:96, j in 1:96)], ct_grid);

# ## The MPI image
#
# 20 × 20 pixels of 2 mm in the MPI device's coordinates. In the MPI device the patient lay shifted
# and rotated by 15 degrees. A registration would compute this relation; here it is given as a
# transform `t`. It maps positions in the CT scanner (the fixed image, whose grid is the output) to
# positions in the MPI device (the moving image, which is resampled), the direction resampling needs
# (see [Grids, Spaces and Transforms](@ref)):

c, s = cos(deg2rad(15)), sin(deg2rad(15))
t = AffineTransform([c -s; s c], [4.0, -3.0]; source=CTScanner(), target=MPIDevice())

# The tracer sits next to the left bone, at `p` in CT coordinates. The MPI image sees it at `t(p)`
# in its own coordinates, blurred to a Gaussian of 3 mm:

p = [-3.0, 4.0]
mpi_grid = GridGeometry((20, 20); spacing=(2.0, 2.0), origin=(-19.0, -19.0), space=MPIDevice())
tracer(x) = exp(-sum(abs2, x - t(p)) / (2 * 3.0^2))
mpi = SpatialImage([tracer(index_to_world(mpi_grid, (i, j))) for i in 1:20, j in 1:20], mpi_grid);

# Both images in the coordinates of their own device. The numbers alone place the tracer at `t(p)`,
# far from the bone it belongs to:

fig = Figure(size=(820, 380)) #hide
lim = ((-25, 25), (-25, 25)) #hide
ax = world_axis(fig, (1, 1); title="CT, in CT scanner coordinates", limits=lim) #hide
world_image!(ax, ct; colorrange=(0, 1)) #hide
ax = world_axis(fig, (1, 2); title="MPI, in MPI device coordinates", limits=lim) #hide
world_image!(ax, mpi; colormap=:inferno, colorrange=(0, 1)) #hide
save("overlay_devices.png", fig) #hide
nothing #hide
# ![The CT and the MPI image, each in its own coordinates](overlay_devices.png)

# ## Resampling through the transform
#
# Without a transform, [`resample`](@ref) refuses, because the two grids are in different spaces:

try #hide
  resample(mpi, ct_grid)
catch e #hide
  print(sprint(showerror, e)) #hide
end #hide

# With the transform, every CT pixel takes the MPI value at its position:

overlay = resample(mpi, ct_grid, t);

# The overlay is on the CT grid, in the CT scanner's space, and its maximum lies at the tracer
# position `p`, up to the 2 mm pixels of the MPI image:

space(overlay), samegrid(geometry(overlay), ct_grid), index_to_world(ct_grid, Tuple(argmax(parent(overlay))))

#-

fig = Figure(size=(440, 400)) #hide
ax = world_axis(fig, (1, 1); title="MPI on the CT grid", limits=lim) #hide
world_image!(ax, ct; colorrange=(0, 1)) #hide
world_image!(ax, overlay; colormap=(:inferno, 0.7), colorrange=(0.1, 1), lowclip=:transparent) #hide
scatter!(ax, [Point2f(p)]; marker=:xcross, color=:cyan, markersize=12) #hide
save("overlay_result.png", fig) #hide
nothing #hide
# ![The MPI image resampled onto the CT grid, over the CT image](overlay_result.png)

# The cross marks `p`. A pixel of the MPI image covers 4 × 4 CT pixels, so the overlay is smooth
# but coarse: resampling places the values, it does not add resolution.
