@testsnippet MPIImage begin
  # (contrast, x, y, z, time) with the grid dimensions in the middle, as MPI reconstructions have
  g = GridGeometry((4, 3, 2); spacing=(0.5, 1.0, 2.0), origin=(1.0, 2.0, 3.0), space=LPS())
  data = reshape(Float32.(1:(2*4*3*2*5)), 2, 4, 3, 2, 5)
  img = TomoImage(data, g; dims=(contrast=CategoricalAxis(["A", "B"]), x=GridAxis(1), y=GridAxis(2),
                                    z=GridAxis(3), time=NumericAxis(0.0:0.25:1.0; unit="s")),
                     meta=(note="test",))
end

@testitem "Axes" tags = [:images] begin
  # numeric and categorical axes are vectors of their values
  t = NumericAxis(0.0:0.5:2.0; unit="s")
  @test t == [0.0, 0.5, 1.0, 1.5, 2.0] && first(t) == 0.0 && TomoImages.unit(t) == "s"
  @test t[2:3] isa NumericAxis && TomoImages.unit(t[2:3]) == "s"
  c = CategoricalAxis(["A", "B"])
  @test c == ["A", "B"] && TomoImages.unit(c) === nothing
  # values: the nearest coordinate on numeric axes, within half a gap at the ends; equal labels
  @test TomoImages.indexof(t, 0.6) == 2
  @test TomoImages.indexof(t, 2.2) == 5
  @test_throws ArgumentError TomoImages.indexof(t, 2.3)
  @test TomoImages.indexof(c, "B") == 2
  # numbers on a categorical axis are labels, not coordinates
  coils = CategoricalAxis([3, 7, 9])
  @test TomoImages.indexof(coils, 7) == 2
  @test_throws ArgumentError TomoImages.indexof(coils, 8)
end

@testitem "TomoImage: construction" tags = [:images] setup = [MPIImage] begin
  @test size(img) == size(data)
  @test parent(img) === data
  @test geometry(img) === g
  @test TomoImages.griddims(img) == (2, 3, 4)
  @test keys(TomoImages.dims(img)) == (:contrast, :x, :y, :z, :time)
  @test space(img) == LPS()
  # without dims: the grid axes x1, x2, ... first, then index vectors
  d = TomoImages.dims(TomoImage(rand(4, 3, 2, 7), g))
  @test keys(d) == (:x1, :x2, :x3, :dim_4) && d.x2 == GridAxis(2) && d.dim_4 == 1:7
  # sizes of data, geometry and axes must agree, and every index axis needs one GridAxis
  @test_throws DimensionMismatch TomoImage(data, g; dims=(x=GridAxis(1), y=GridAxis(2), z=GridAxis(3),
                                                             a=1:2, b=1:5))
  @test_throws ArgumentError TomoImage(data, g; dims=(a=1:2, x=GridAxis(1), y=GridAxis(1), z=GridAxis(3), b=1:5))
  @test_throws DimensionMismatch TomoImage(data, g; dims=(contrast=["A"], x=GridAxis(1), y=GridAxis(2),
                                                             z=GridAxis(3), time=1:5))
  @test_throws ArgumentError TomoImage(data, g; dims=(contrast=1:2, x=GridAxis(1)))
end

@testitem "TomoImage: selections keep the geometry" tags = [:images] setup = [MPIImage] begin
  # values are those of the parent; voxel positions are those of the full grid
  s = img[:, 2:4, :, 2, 1:2:5]
  @test parent(s) == data[:, 2:4, :, 2, 1:2:5]
  @test TomoImages.griddims(s) == (2, 3)
  @test index_to_world(geometry(s), (1, 1)) == index_to_world(g, (2, 1, 2))
  @test TomoImages.dims(s).time == [0.0, 0.5, 1.0]
  # an integer drops the dimension; a dropped grid axis renumbers the remaining ones
  v = view(img, 1, :, 2, :, 3)
  @test v isa TomoImage && keys(TomoImages.dims(v)) == (:x, :z)
  @test TomoImages.dims(v).z == GridAxis(2) && parent(v) == data[1, :, 2, :, 3]
  # by name, integers and ranges are indices
  @test parent(img[contrast=2, time=3]) == data[2, :, :, :, 3]
  @test parent(img[x=2:3]) == data[:, 2:3, :, :, :]
  @test_throws ArgumentError img[time=0.5]
  @test_throws ArgumentError img[echo=1]
  # by value with selectvalue: the nearest time within half a step, the equal label
  @test parent(selectvalue(img; contrast="B", time=0.5)) == data[2, :, :, :, 3]
  @test parent(selectvalue(img; contrast="A", time=0.6)) == data[1, :, :, :, 3]
  @test_throws ArgumentError selectvalue(img; time=1.2)
  @test_throws ArgumentError selectvalue(img; contrast="C")
  @test_throws ArgumentError selectvalue(img; x=1.0)
  # integers on all grid dimensions leave a plain array
  @test img[1, 1, 1, 1, :] == data[1, 1, 1, 1, :]
  @test img[1, 2, 3, 2, 4] == data[1, 2, 3, 2, 4]
end

@testitem "TomoImage: permutedims and reverse" tags = [:images] setup = [MPIImage] begin
  # permutedims reorders data and dims; the geometry and every voxel position stay
  p = permutedims(img, (5, 4, 3, 2, 1))
  @test parent(p) == permutedims(data, (5, 4, 3, 2, 1))
  @test geometry(p) === g && TomoImages.griddims(p) == (4, 3, 2)
  @test keys(TomoImages.dims(p)) == (:time, :z, :y, :x, :contrast)
  # reversing a grid dimension keeps every voxel at its position: the geometry flips with it
  r = reverse(img; dims=2)
  @test parent(r) == reverse(data; dims=2)
  for i in 1:4
    @test index_to_world(geometry(r), (i, 1, 1)) == index_to_world(g, (5 - i, 1, 1))
  end
  # reversing a value axis reverses its values
  @test TomoImages.dims(reverse(img; dims=5)).time == [1.0, 0.75, 0.5, 0.25, 0.0]
end

@testitem "TomoImage: broadcasting and similar" tags = [:images] setup = [MPIImage] begin
  b = 2 .* img .+ 1
  @test b isa TomoImage && geometry(b) === g && TomoImages.dims(b) === TomoImages.dims(img) && meta(b) === meta(img)
  @test parent(b) == 2 .* data .+ 1
  c = similar(img)
  c .= img .^ 2
  @test parent(c) == data .^ 2 && geometry(c) === g
  @test copy(img) == img && parent(copy(img)) !== data
end

@testitem "TomoImage: broadcasting checks space and grid" tags = [:images] begin
  g = GridGeometry((4, 3); spacing=(1.0, 1.0), space=LPS())
  a = TomoImage(rand(4, 3), g)
  # images on the same grid combine, the result keeps the grid
  @test geometry(a .+ TomoImage(rand(4, 3), g)) === g
  # another space or another grid of the same size is an error, not a silent choice of one geometry
  @test_throws ArgumentError a .+ TomoImage(rand(4, 3), GridGeometry((4, 3); space=RAS()))
  @test_throws ArgumentError a .+ TomoImage(rand(4, 3), GridGeometry((4, 3); origin=(0.5, 0.0), space=LPS()))
  b = similar(a)
  @test_throws ArgumentError b .= TomoImage(rand(4, 3), GridGeometry((4, 3); space=RAS()))
end
