using Test
using Integration

linear(x) = 2x - 1
intlinear(a, b) = b*(b - 1) - a * (a - 1)
parabolic(x) = 3x^2 - 2x + 1
intparabolic(a, b) = (b^3 - b^2 + b) - (a^3 - a^2 + a)
expsin(x) = x * exp(sin(2x))

atanintegrand(x) = 1 / (1 + x^2)
intatan() = 2 * atan(5)

@testset "Integration" begin

@testset "Midpoint" begin
    a, b = -2, 1
    @testset "Квадратурная формула" begin
        @test integrate(linear, a, b; method=Midpoint()) == intlinear(a, b)
        @test integrate(parabolic, a, b; method=Midpoint()) == 3 * parabolic(-0.5)
        @test integrate(expsin, a, b; method=Midpoint()) == 3 * expsin(-0.5)
        let
            Ia0 = integrate(linear, a, 0; method=Midpoint())
            I0b = integrate(linear, 0, b; method=Midpoint())
            @test Ia0 + I0b == intlinear(a, b)
        end
    end
    @testset "Составная квадратурная формула" begin
        @test integrate(linear, a, b, 10; method=Midpoint()) == intlinear(a, b)
        @test integrate(atanintegrand, -5, 5, 1000; method=Midpoint()) ≈ intatan() atol=1e-6
    end
end

@testset "Trapezoid" begin
    a, b = -2, 1
    @testset "Квадратурная формула" begin
        @test integrate(linear, a, b; method=Trapezoid()) == intlinear(a, b)
        @test integrate(parabolic, a, b; method=Trapezoid()) == 3 * sum(parabolic, (a, b)) / 2
        @test integrate(expsin, a, b; method=Trapezoid()) == 3 * sum(expsin, (a, b)) / 2
        let
            Ia0 = integrate(linear, a, 0; method=Trapezoid())
            I0b = integrate(linear, 0, b; method=Trapezoid())
            @test Ia0 + I0b == intlinear(a, b)
        end
    end
    @testset "Составная квадратурная формула" begin
        @test integrate(linear, a, b, 10; method=Trapezoid()) == intlinear(a, b)
        @test integrate(atanintegrand, -5, 5, 1000; method=Trapezoid()) ≈ intatan() atol=1e-6
    end
end

@testset "Simpson" begin
    a, b = -2, 1
    @testset "Квадратурная формула" begin
        @test integrate(linear, a, b; method=Simpson()) == intlinear(a, b)
        @test integrate(parabolic, a, b; method=Simpson()) == intparabolic(a, b)
    end
    @testset "Составная квадратурная формула" begin
        @test integrate(parabolic, a, b, 10; method=Simpson()) == intparabolic(a, b)
        @test integrate(atanintegrand, -5, 5, 100; method=Simpson()) ≈ intatan() atol=1e-6
    end
end

@testset "Gauss" begin
    @testset "Интегрирование полинома" begin
        @test integrate(x->14*x^13, -1, 2; method=Gauss()) ≈ 2^14 - 1
        @test integrate(x->15*x^14, -1, 2; method=Gauss()) ≈ 2^15 - 1 atol=1
    end
    @testset "Интегрирование 1/(1+x²)" begin
        a, b = -5, 5
        @test integrate(atanintegrand, a, b; method=Gauss()) ≈ 3.0806104010709623 atol=1e-6
        @testset "Составная квадратурная формула" begin
            @test integrate(atanintegrand, a, b, 6; method=Gauss()) ≈ intatan() rtol=1e-6
        end
    end
end

@testset "Kronrod" begin
    @testset "Интегрирование полинома" begin
        @test integrate(x->14*x^13, -1, 2; method=Kronrod()) ≈ 2^14 - 1
        @test integrate(x->15*x^14, -1, 2; method=Kronrod()) ≈ 2^15 - 1 atol=3
        @test integrate(x->30*x^29, -1, 2; method=Kronrod()) ≈ 2^30 - 1 rtol=1e-6
    end
    @testset "Интегрирование 1/(1+x²)" begin
        a, b = -5, 5
        @test integrate(atanintegrand, a, b; method=Kronrod()) ≈ 2.7631456512762456 atol=1e-6
        @testset "Составная квадратурная формула" begin
            @test integrate(atanintegrand, a, b, 4; method=Kronrod()) ≈ intatan() rtol=1e-6
        end
    end
end

end # testset Integration
