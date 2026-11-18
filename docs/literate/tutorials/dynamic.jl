# # Analyze a Dynamic Image
#
# This tutorial is for users of reconstructed images. A dynamic image shows a tracer passing
# through two vessels over time; a label map drawn on a finer anatomical grid marks the vessels.
# The tutorial selects frames by time, brings the label map onto the grid of the dynamic image and
# computes the mean concentration in each vessel over time.

using TomoImages, CairoMakie
include(joinpath(pkgdir(TomoImages), "docs", "utils", "figures.jl")) #hide
nothing #hide

# ## The images
#
# The dynamic image: 32 × 32 pixels of 2 mm and 60 frames of 0.05 s. The tracer passes vessel A
# (a disk of radius 6 mm) first and vessel B (radius 5 mm) later:

g = GridGeometry((32, 32); spacing=(2.0, 2.0), origin=(-31.0, -31.0))
times = NumericAxis(0.05 .* (0:59); unit="s")
bolus(t, t0) = t > t0 ? ((t - t0) / 0.4)^2 * exp(-(t - t0) / 0.4) : 0.0
vessel_a(x, y) = hypot(x + 10, y) < 6
vessel_b(x, y) = hypot(x - 12, y - 4) < 5
function concentration(x, y, t)
  vessel_a(x, y) * bolus(t, 0.2) + vessel_b(x, y) * 0.6 * bolus(t, 0.8)
end
data = [concentration(index_to_world(g, (i, j))..., t) for i in 1:32, j in 1:32, t in times]
dyn = SpatialImage(data, g; dims=(x=GridAxis(1), y=GridAxis(2), time=times));

# The label map was drawn on an anatomical image of 0.5 mm pixels, in the same space: label 1 for
# vessel A, label 2 for vessel B, each drawn 2 mm inside the vessel wall:

fine = GridGeometry((128, 128); spacing=(0.5, 0.5), origin=(-31.75, -31.75))
label(x, y) = hypot(x + 10, y) < 4 ? 0x01 : hypot(x - 12, y - 4) < 3 ? 0x02 : 0x00
labels = SpatialImage([label(index_to_world(fine, (i, j))...) for i in 1:128, j in 1:128], fine);

# ## Frames by time
#
# [`selectvalue`](@ref) picks the frame nearest to a time; the result is a 2D image on the same
# grid. Integers select by index: `dyn[time=21]` is the same frame as time 1.0 s:

frame = selectvalue(dyn; time=1.0)
parent(frame) == parent(dyn[time=21])

#-

fig = Figure(size=(820, 380)) #hide
lim = ((-32, 32), (-32, 32)) #hide
ax = world_axis(fig, (1, 1); title="frame at 1.0 s", limits=lim) #hide
world_image!(ax, frame; colormap=:inferno, colorrange=(0, 0.6)) #hide
ax = world_axis(fig, (1, 2); title="label map (0.5 mm)", limits=lim) #hide
world_image!(ax, labels; colormap=[:black, :steelblue, :darkorange], colorrange=(0, 2)) #hide
save("dynamic_images.png", fig) #hide
nothing #hide
# ![A frame of the dynamic image and the label map](dynamic_images.png)

# ## The label map on the dynamic grid
#
# The two images are in the same space but on different grids, so they cannot be combined pixel by
# pixel:

isapprox(geometry(labels), geometry(dyn))

# [`resample`](@ref) brings the label map onto the grid of the dynamic image. Labels are integers,
# so the default resampler is [`NearestResampler`](@ref): every 2 mm pixel takes the label at its
# center, and no label is averaged with its neighbors:

on_grid = resample(labels, geometry(dyn))
eltype(on_grid), sort(unique(parent(on_grid)))

# ## Concentration over time
#
# For every label, the mean of the dynamic image over its pixels, frame by frame:

mean_in(img, mask) = [sum(parent(img)[:, :, k][mask]) / count(mask) for k in 1:size(img, 3)]
curve_a = mean_in(dyn, parent(on_grid) .== 1)
curve_b = mean_in(dyn, parent(on_grid) .== 2);

# The labels lie inside the vessels, so their pixels hold exactly the simulated concentrations, and
# the curves reproduce them up to the rounding of the mean:

maximum(abs, curve_a .- bolus.(times, 0.2)), maximum(abs, curve_b .- 0.6 .* bolus.(times, 0.8))

#-

fig = Figure(size=(620, 320)) #hide
ax = Axis(fig[1, 1]; xlabel="time (s)", ylabel="mean concentration") #hide
lines!(ax, collect(times), curve_a; label="vessel A") #hide
lines!(ax, collect(times), curve_b; label="vessel B") #hide
axislegend(ax) #hide
save("dynamic_curves.png", fig) #hide
nothing #hide
# ![Mean concentration in the two vessels over time](dynamic_curves.png)

# The time axis of the image gives the times of the frames, so the curves are plotted against
# seconds without separate bookkeeping.
