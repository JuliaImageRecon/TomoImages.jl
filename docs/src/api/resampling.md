# Resampling

[`resample`](@ref) puts an image onto another grid, in the same space or through a transform. The
resampler is its first argument; [Write a Resampler](@ref) shows how to add one.

```@docs
resample
AbstractResampler
NearestResampler
LinearResampler
BSplineResampler
TomoImages.default_resampler
TomoImages.supports
TomoImages.resample!
TomoImages.resample_eltype
TomoImages.IndexMap
TomoImages.index_map
```
