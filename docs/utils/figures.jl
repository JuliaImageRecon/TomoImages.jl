# Figures of the documentation (CairoMakie): 2D images drawn in the coordinates of their space,
# every pixel as the quadrilateral it covers, so that rotated grids appear rotated. Included
# invisibly by the pages.

using CairoMakie, TomoImages

# the four corners of pixel (i, j) in the space of g, as the world coordinates `coords` (two of them)
pixel_corners(g, i, j, coords) = [Point2f(index_to_world(g, (i + a, j + b))[collect(coords)])
                                  for (a, b) in ((-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5))]

# a 2D image (or a slice of a 3D image) in the world coordinates `coords`
function world_image!(ax, img; coords=(1, 2), colormap=:grays, colorrange=extrema(parent(img)), kwargs...)
  g = geometry(img)
  D = TomoImages.direction(g)
  if size(D, 1) == 2 && D == [1 0; 0 1]
    return heatmap!(ax, axiscoords(g, 1), axiscoords(g, 2), parent(img); colormap, colorrange, kwargs...)
  end
  quads = [pixel_corners(g, i, j, coords) for i in 1:size(g, 1), j in 1:size(g, 2)]
  poly!(ax, vec(quads); color=vec(Float32.(parent(img))), colormap, colorrange, strokewidth=0, kwargs...)
end

# the outline of a 2D grid: the outer edges of its pixels
function outline!(ax, g; kwargs...)
  n1, n2 = size(g)
  corners = [index_to_world(g, c) for c in ((0.5, 0.5), (n1 + 0.5, 0.5), (n1 + 0.5, n2 + 0.5), (0.5, n2 + 0.5), (0.5, 0.5))]
  lines!(ax, Point2f.(corners); kwargs...)
end

function world_axis(fig, pos; title="", limits=nothing, xlabel="x (mm)", ylabel="y (mm)")
  ax = Axis(fig[pos...]; aspect=DataAspect(), title, xlabel, ylabel, xgridvisible=false, ygridvisible=false)
  limits === nothing || limits!(ax, limits...)
  ax
end

# a grid with its pixel edges, pixel centers, origin and direction arrows
function grid_figure(g)
  fig = Figure(size=(560, 460))
  ax = world_axis(fig, (1, 1); title="grid in its space")
  n1, n2 = size(g)
  for i in 0.5:1:n1+0.5
    lines!(ax, Point2f.([index_to_world(g, (i, 0.5)), index_to_world(g, (i, n2 + 0.5))]); color=:gray80)
  end
  for j in 0.5:1:n2+0.5
    lines!(ax, Point2f.([index_to_world(g, (0.5, j)), index_to_world(g, (n1 + 0.5, j))]); color=:gray80)
  end
  centers = [Point2f(index_to_world(g, (i, j))) for i in 1:n1, j in 1:n2]
  scatter!(ax, vec(centers); color=:gray50, markersize=5)
  o = Point2f(TomoImages.origin(g))
  scatter!(ax, [o]; color=:darkorange, markersize=12)
  text!(ax, o; text="origin: center of pixel (1, 1)", offset=(8, -16), color=:darkorange, fontsize=13)
  for (d, c) in ((1, :steelblue), (2, :seagreen))
    v = Point2f(TomoImages.direction(g)[:, d] * spacing(g)[d] * 2)
    arrows2d!(ax, [o], [v]; color=c)
    text!(ax, o + v; text="index axis $d", color=c, fontsize=13, offset=(4, 4))
  end
  fig
end

# the outline of a mask on an axis-aligned 2D grid, in world coordinates
mask_outline!(ax, mask; kwargs...) = contour!(ax, axiscoords(geometry(mask), 1), axiscoords(geometry(mask), 2),
                                              Float32.(parent(mask)); levels=[0.5], kwargs...)
