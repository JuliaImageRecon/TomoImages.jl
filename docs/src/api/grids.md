# Spaces and Grids

A grid geometry places voxels in a coordinate space: voxel `i` has its center at
`origin + direction * (spacing .* (i .- 1))`. The space names the coordinate system these positions
are written in. [Grids, Spaces and Transforms](@ref) explains the conventions.

```@docs
TomoImages
```

## Spaces

The coordinate systems that geometries and transforms refer to.

```@docs
AbstractSpace
DefaultSpace
LPS
RAS
TomoImages.spacename
```

## Grid geometries

A grid geometry, its parts, and the functions that map between voxel indices and positions.
`origin` and `direction` are public but not exported, since the names are common.

```@docs
AbstractGeometry
GridGeometry
spacing
TomoImages.origin
TomoImages.direction
space
TomoImages.worlddims
index_to_world
world_to_index
affine
fov
fov_center
axiscoords
Base.isapprox(::GridGeometry, ::GridGeometry)
Base.in(::Any, ::GridGeometry)
TomoImages.subgrid
```
