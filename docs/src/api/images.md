# Images

An image is an array whose grid dimensions lie on a grid geometry, with a named description of
every dimension.

```@docs
AbstractSpatialImage
SpatialImage
geometry
TomoImages.dims
TomoImages.spatialdims
meta
selectvalue
```

## Axes

Every dimension of an image is a grid axis or is described by a vector of values.

```@docs
GridAxis
NumericAxis
CategoricalAxis
TomoImages.unit
TomoImages.indexof
```
