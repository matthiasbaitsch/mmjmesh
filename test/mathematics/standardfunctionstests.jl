@testitem "Standard functions" setup = [Validate] begin

    using Test

    using MMJMesh
    using MMJMesh.Mathematics


    # # Sine and Cosine

    s = Sin(0 .. 2π)
    c = Cos(0 .. 2π)
    @test Validate.validate(s)
    @test Validate.validate(antiderivative(s))
    @test Validate.validate(c)
    @test Validate.validate(antiderivative(c))
    @test antiderivative(Sin(), 2) == -Sin()

    x = parameter(0 .. 2)
    @test sin(x) == Sin(0 .. 2)
    @test cos(x) == Cos(0 .. 2)


    # # Exp

    f = Exp(0 .. 2π)
    @test Validate.validate(f)
    @test Validate.validate(antiderivative(f))
    @test exp(x) == Exp(0 .. 2)


    # # Sinh and Cosh

    sh = Sinh(0 .. 2)
    ch = Cosh(0 .. 2)
    @test Validate.validate(sh)
    @test Validate.validate(antiderivative(sh))
    @test Validate.validate(ch)
    @test Validate.validate(antiderivative(ch))
    @test antiderivative(Sinh(), 2) == Sinh()
    @test sinh(x) == Sinh(0 .. 2)
    @test cosh(x) == Cosh(0 .. 2)


    # # Abs

    f = Abs(-1 .. 2)
    @test f(-1) == 1
    @test f(2) == 2
    @test derivativeat(f, -0.5) == -1
    @test derivativeat(f, 0.5) == 1
    @test derivativeat(f, 0.5, 2) == 0
    @test antiderivative(f)(-1) ≈ -0.5
    @test antiderivative(f)(2) ≈ 2
    @test integrate(f, -1 .. 2) ≈ 2.5
    @test abs(x) == Abs(0 .. 2)
    @test pois(f) == [0]
    @test pois(Abs(0.5 .. 2)) == []
    @test Abs()' == Sign()
    @test Validate.validate(Abs(0.5 .. 2), atol=1e-6)
    @test Validate.validate(Abs(-2 .. -0.5), atol=1e-6)


    # # Sign

    f = Sign(-1 .. 2)
    @test f(-1) == -1
    @test f(2) == 1
    @test derivativeat(f, 0.5) == 0
    @test antiderivative(f) == Abs(-1 .. 2)
    @test sign(x) == Sign(0 .. 2)
    @test pois(f) == [0]
    @test pois(Sign(0.5 .. 2)) == []
    @test Validate.validate(Sign(0.5 .. 2), atol=1e-6)
    @test Validate.validate(Sign(-2 .. -0.5), atol=1e-6)
end
