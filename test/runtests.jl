using Test
using Integration

@testset "Midpoint" begin
    integrate(x -> x - 1, -2, 1; method=Midpoint())
end
