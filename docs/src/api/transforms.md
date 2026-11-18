# Transforms

Transforms map positions from one space to another and carry both spaces, so that composition and
resampling check them.

```@docs
AbstractSpatialTransform
AffineTransform
source
target
compose
transform_between
transform_geometry
TomoImages.checkspaces
```
