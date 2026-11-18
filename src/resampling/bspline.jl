supports(::BSplineResampler, img) = parent(img) isa Array

bspline(degree) = degree == 0 ? Interpolations.BSpline(Interpolations.Constant()) :
                  degree == 1 ? Interpolations.BSpline(Interpolations.Linear()) :
                  degree == 2 ? Interpolations.BSpline(Interpolations.Quadratic(Interpolations.Line(Interpolations.OnGrid()))) :
                  Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid())))

# One interpolant per non-spatial index; dimensions of length 1 are not interpolated.
function resample!(r::BSplineResampler, out::AbstractArray, img::AbstractSpatialImage, imap::IndexMap;
                   fill=zero(eltype(out)))
  data = parent(img)
  sd = spatialdims(img)
  others = Tuple(d for d in 1:ndims(data) if !(d in sd))
  n = map(d -> size(data, d), sd)
  kinds = map(m -> m == 1 ? Interpolations.NoInterp() : bspline(r.degree), n)
  for e in CartesianIndices(map(d -> size(data, d), others))
    I = ntuple(d -> d in sd ? Colon() : e[findfirst(==(d), others)], ndims(data))
    itp = Interpolations.interpolate(view(data, I...), kinds)
    target = view(out, I...)
    for j in CartesianIndices(map(d -> size(out, d), sd))
      p = imap(Tuple(j))
      target[j] = inside(p, n) ? itp(ntuple(k -> n[k] == 1 ? 1 : clamp(p[k], 1, n[k]), length(n))...) : fill
    end
  end
  out
end
