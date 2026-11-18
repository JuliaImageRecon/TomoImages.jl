# Spaces and Grids

A grid geometry places voxels in a coordinate space: voxel `i` has its center at
`origin + direction * (spacing .* (i .- 1))`. The space names the coordinate system these positions
are written in. [Grids, Spaces and Transforms](@ref) explains the conventions.

```@docs
TomoImages
```

## Spaces

```@docs
AbstractSpace
DefaultSpace
LPS
RAS
spacename
```

## Grid geometries

```@docs
AbstractGeometry
GridGeometry
spacing
origin
direction
space
TomoImages.worlddims
index_to_world
world_to_index
affine
fov
fov_center
axiscoords
samegrid
subgrid
```
