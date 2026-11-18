# # Explore a Spatial Image
#
# This tutorial takes one 3D image and reads what it says about itself: where its voxels are, how
# to show it in millimeters, how to find the voxel at a position and the position of a voxel, how to
# measure in it, and how to look at it along a tilted structure.

using TomoImages, CairoMakie
include(joinpath(pkgdir(TomoImages), "docs", "utils", "figures.jl")) #hide
nothing #hide

# ## The image
#
# A volume of 64 × 64 × 40 voxels in the patient space [`LPS`](@ref), with voxels of 1 × 1 × 1.5 mm.
# It shows a body (value 0.2), two balls (1.0 and 0.8) and a thin rod (1.2) that runs obliquely
# through the volume. The values come from a function of the position in millimeters:

g = GridGeometry((64, 64, 40); spacing=(1.0, 1.0, 1.5), origin=(-31.5, -31.5, -29.25), space=LPS())
rod_axis = [1.0, 0.5, 1.0] / 1.5          # a unit vector
function phantom(x, y, z)
  p = [x, y, z]
  t = p' * rod_axis                       # distance along the rod from its center at the origin
  abs(t) < 25 && sum(abs2, p - t * rod_axis) < 2^2 && return 1.2
  sum(abs2, p - [-15, 10, -10]) < 5^2 && return 1.0
  sum(abs2, p - [12, -8, 12]) < 4^2 && return 0.8
  (x / 28)^2 + (y / 25)^2 + (z / 28)^2 < 1 && return 0.2
  0.0
end
img = SpatialImage([phantom(index_to_world(g, (i, j, k))...) for i in 1:64, j in 1:64, k in 1:40], g;
                   dims=(x=GridAxis(1), y=GridAxis(2), z=GridAxis(3)));

# The image knows its size, its geometry and the names of its dimensions:

summary(img)

#-

geometry(img)

# The origin is the center of the first voxel; the spacing is the distance between voxel centers;
# the direction says that the index axes run along the axes of the space. The volume extends from
# edge to edge over [`fov`](@ref) millimeters around its center:

fov(g), fov_center(g)

# ## Slices in millimeters
#
# Selecting one index of a dimension gives a 2D image that keeps its place in the volume. The three
# slices through the center, drawn in the coordinates of the space, so that the voxels of 1.5 mm
# along `z` appear as long as they are:

axial, coronal, sagittal = img[z=20], img[y=32], img[x=32]
index_to_world(geometry(axial), (1, 1))

#-

fig = Figure(size=(1050, 380)) #hide
for (j, (title, s, coords, xl, yl)) in enumerate((("axial (z = 20)", axial, (1, 2), "x (mm)", "y (mm)"), #hide
                                                 ("coronal (y = 32)", coronal, (1, 3), "x (mm)", "z (mm)"), #hide
                                                 ("sagittal (x = 32)", sagittal, (2, 3), "y (mm)", "z (mm)"))) #hide
  ax = world_axis(fig, (1, j); title, xlabel=xl, ylabel=yl, limits=((-33, 33), (-33, 33))) #hide
  world_image!(ax, s; coords, colorrange=(0, 1.2)) #hide
end #hide
save("explore_slices.png", fig) #hide
nothing #hide
# ![Axial, coronal and sagittal slices through the center, in millimeters](explore_slices.png)

# ## Positions and voxels
#
# [`index_to_world`](@ref) gives the position of a voxel, [`world_to_index`](@ref) the voxel at a
# position. The center of the second ball lies at `(12, -8, 12)` mm; the voxel there and its value:

i = world_to_index(g, (12.0, -8.0, 12.0))

#-

img[round.(Int, i)...]

# Positions outside the voxels of the image are not in the grid:

(12.0, -8.0, 12.0) in g, (0.0, 0.0, 40.0) in g

# ## Measuring in millimeters
#
# The first ball, as the set of its voxels: their positions give its center and its volume in the
# units of the space, independent of the voxel size along each axis:

ball = findall(==(1.0), parent(img))
positions = [index_to_world(g, Tuple(c)) for c in ball]
center = sum(positions) / length(positions)
volume = length(ball) * prod(spacing(g))
center, volume, 4 / 3 * pi * 5^3

# The center lies within a fraction of a voxel of `(-15, 10, -10)` mm, where the ball was placed;
# the volume of its voxels is that of a ball of radius 5 mm, up to the voxels at its surface.
#
# ## Along the rod
#
# The rod runs obliquely through all three slices above, so no slice shows it whole. A grid whose
# first index axis runs along the rod shows it: 80 × 40 × 1 voxels of 0.75 mm in the plane spanned
# by the rod and the `z`-perpendicular direction, centered at the origin. Its direction matrix holds
# the rod axis, a perpendicular direction and their cross product:

u = rod_axis
v = [-u[2], u[1], 0.0] / hypot(u[1], u[2])
w = [u[2] * v[3] - u[3] * v[2], u[3] * v[1] - u[1] * v[3], u[1] * v[2] - u[2] * v[1]]
D = hcat(u, v, w)
plane = GridGeometry((80, 40, 1); spacing=(0.75, 0.75, 0.75), origin=-D * [0.75 * 79 / 2, 0.75 * 39 / 2, 0.0],
                     direction=D, space=LPS())
reformat = resample(img, plane);

# A profile along the rod is the same with a grid of one voxel across: 120 points of 0.5 mm along
# the rod axis:

line = GridGeometry((120, 1, 1); spacing=(0.5, 0.5, 0.5), origin=-u * 0.5 * 119 / 2, direction=D, space=LPS())
profile = parent(resample(img, line))[:, 1, 1];

#-

fig = Figure(size=(1000, 360)) #hide
ax = Axis(fig[1, 1]; aspect=DataAspect(), title="plane along the rod", xlabel="along the rod (mm)", #hide
          ylabel="across (mm)") #hide
heatmap!(ax, 0.75 .* ((1:80) .- 40.5), 0.75 .* ((1:40) .- 20.5), parent(reformat)[:, :, 1]; colormap=:grays, #hide
         colorrange=(0, 1.2)) #hide
ax = Axis(fig[1, 2]; title="profile along the rod", xlabel="along the rod (mm)", ylabel="value") #hide
lines!(ax, 0.5 .* ((1:120) .- 60.5), profile) #hide
save("explore_rod.png", fig) #hide
nothing #hide
# ![The plane along the rod and a profile along its axis](explore_rod.png)

# In the plane the rod appears whole, 50 mm long; the profile shows its value 1.2 over that length
# and the body around it. Both are resampled with multilinear interpolation, so the edges of the
# rod are smoothed over about one voxel.
