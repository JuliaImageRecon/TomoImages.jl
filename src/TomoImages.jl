"""
    TomoImages

Images that know where they are in space: grid geometries with an origin, a spacing, a direction and
a named coordinate space; images on such grids with named non-spatial axes; transforms between
spaces; and resampling from one grid onto another, on CPUs and GPUs.
"""
module TomoImages

using LinearAlgebra, StaticArrays
using Adapt
using KernelAbstractions
import Interpolations

const KA = KernelAbstractions

export AbstractSpace, DefaultSpace, LPS, RAS, spacename
export AbstractGeometry, GridGeometry, spacing, origin, direction, space, affine, index_to_world,
       world_to_index, fov, fov_center, axiscoords, samegrid, subgrid
export AbstractSpatialTransform, AffineTransform, source, target, compose, transform_between, transform_geometry
export AbstractAxis, NumericAxis, CategoricalAxis
export AbstractSpatialImage, SpatialImage, geometry, spatialdims, extra, meta
export AbstractResampler, NearestResampler, LinearResampler, BSplineResampler, resample, resample!, default_resampler, supports

public IndexMap, index_map, resample_eltype, checkspaces, worlddims, coords, labels, unit, indexof

include("spaces/spaces.jl")
include("geometry/geometry.jl")
include("transforms/transforms.jl")
include("images/images.jl")
include("resampling/resampling.jl")

end
