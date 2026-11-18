supports(::BSplineResampler, img) = parent(img) isa Array

# linear extrapolation as the boundary condition of the prefilter (degree 2 and 3)
const ONGRID = Interpolations.Line(Interpolations.OnGrid())
bspline(degree) = degree == 0 ? Interpolations.BSpline(Interpolations.Constant()) :
                  degree == 1 ? Interpolations.BSpline(Interpolations.Linear()) :
                  degree == 2 ? Interpolations.BSpline(Interpolations.Quadratic(ONGRID)) :
                  Interpolations.BSpline(Interpolations.Cubic(ONGRID))

# One interpolant per non-grid index; dimensions of length 1 are not interpolated. The views have
# the grid dimensions in array order, the index map in the order of the geometry's index axes:
# view dimension m is index axis order[m].
function resample!(r::BSplineResampler, out::AbstractArray, img::AbstractTomoImage, imap::IndexMap;
                   fillvalue=zero(eltype(out)))
  data = parent(img)
  sd = griddims(img)
  order = Tuple(sortperm(collect(sd)))
  others = Tuple(d for d in 1:ndims(data) if !(d in sd))
  n = map(d -> size(data, d), sd)
  kinds = map(k -> n[k] == 1 ? Interpolations.NoInterp() : bspline(r.degree), order)
  for e in CartesianIndices(map(d -> size(data, d), others))
    I = ntuple(d -> d in sd ? Colon() : e[findfirst(==(d), others)], ndims(data))
    itp = Interpolations.interpolate(view(data, I...), kinds)
    target = view(out, I...)
    for j in CartesianIndices(map(d -> size(out, d), sd))
      p = imap(Tuple(j))
      inside(p, n) || (target[ntuple(m -> j[order[m]], length(order))...] = fillvalue; continue)
      q = ntuple(m -> (k = order[m]; n[k] == 1 ? 1 : clamp(p[k], 1, n[k])), length(order))
      target[ntuple(m -> j[order[m]], length(order))...] = itp(q...)
    end
  end
  out
end
