@testsnippet MPIImage begin
  # (contrast, x, y, z, time) with the spatial dimensions in the middle, as MPI reconstructions have
  g = GridGeometry((4, 3, 2); spacing=(0.5, 1.0, 2.0), origin=(1.0, 2.0, 3.0), space=LPS())
  data = reshape(Float32.(1:(2*4*3*2*5)), 2, 4, 3, 2, 5)
  img = SpatialImage(data, g; spatialdims=(2, 3, 4),
                     extra=(contrast=CategoricalAxis(["A", "B"]), time=NumericAxis(0.0:0.25:1.0; unit="s")),
                     meta=(note="test",))
end

@testitem "SpatialImage: construction" tags = [:images] setup = [MPIImage] begin
  @test size(img) == size(data)
  @test parent(img) === data
  @test geometry(img) === g
  @test spatialdims(img) == (2, 3, 4)
  @test keys(extra(img)) == (:contrast, :time)
  @test space(img) == LPS()
  # sizes of data, geometry and axes must agree
  @test_throws DimensionMismatch SpatialImage(data, g; spatialdims=(1, 2, 3),
                                              extra=(a=CategoricalAxis(1:2), b=NumericAxis(1:5)))
  @test_throws ArgumentError SpatialImage(data, g; spatialdims=(2, 3, 4))
  @test_throws DimensionMismatch SpatialImage(data, g; spatialdims=(2, 3, 4),
                                              extra=(contrast=CategoricalAxis(["A"]), time=NumericAxis(1:5)))
end

@testitem "SpatialImage: selections keep the geometry" tags = [:images] setup = [MPIImage] begin
  using StaticArrays
  # values are those of the parent; voxel positions are those of the full grid
  s = img[:, 2:4, :, 2, 1:2:5]
  @test parent(s) == data[:, 2:4, :, 2, 1:2:5]
  @test spatialdims(s) == (2, 3)
  @test index_to_world(geometry(s), (1, 1)) == index_to_world(g, (2, 1, 2))
  @test TomoImages.coords(extra(s).time) == [0.0, 0.5, 1.0]
  # an integer on an extra axis drops the axis
  v = view(img, 1, :, :, :, 3)
  @test v isa SpatialImage && spatialdims(v) == (1, 2, 3) && keys(extra(v)) == ()
  @test parent(v) == data[1, :, :, :, 3]
  # selecting by value on the extra axes
  b = img[contrast="B", time=0.5]
  @test parent(b) == data[2, :, :, :, 3]
  @test geometry(b) === g
  # the nearest time within half a step; beyond the ends it is an error
  @test parent(img[contrast="A", time=0.6]) == data[1, :, :, :, 3]
  @test_throws ArgumentError img[time=1.2]
  @test_throws ArgumentError img[contrast="C"]
  @test_throws ArgumentError img[echo=1]
  # integers on all spatial dimensions leave a plain array
  @test img[1, 1, 1, 1, :] == data[1, 1, 1, 1, :]
  @test img[1, 2, 3, 2, 4] == data[1, 2, 3, 2, 4]
end

@testitem "SpatialImage: broadcasting and similar" tags = [:images] setup = [MPIImage] begin
  b = 2 .* img .+ 1
  @test b isa SpatialImage && geometry(b) === g && extra(b) === extra(img) && meta(b) === meta(img)
  @test parent(b) == 2 .* data .+ 1
  c = similar(img)
  c .= img .^ 2
  @test parent(c) == data .^ 2 && geometry(c) === g
  @test copy(img) == img && parent(copy(img)) !== data
end
