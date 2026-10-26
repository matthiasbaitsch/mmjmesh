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
end
