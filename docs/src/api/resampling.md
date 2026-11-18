# Resampling

[`resample`](@ref) puts an image onto another grid, in the same space or through a transform. The
resampler is its first argument.

```@docs
resample
NearestResampler
LinearResampler
BSplineResampler
TomoImages.default_resampler
```

## Implementing a resampler

A resampler is a type of its own with a method of `TomoImages.resample!`; [Write a Resampler](@ref)
shows one.

```@docs
AbstractResampler
TomoImages.resample!
TomoImages.supports
TomoImages.resample_eltype
TomoImages.IndexMap
TomoImages.index_map
```
