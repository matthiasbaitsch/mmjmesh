@testitem "Symbolics" begin

    using Test
    using Symbolics

    using MMJMesh
    using MMJMesh.MMJBase

    @variables a, b
    u = [2.0 + 3.0a, b]
    v = [0.5a, 0.348b]
    w = [0.5a 0.348b; 3b a/2]

    # # Make rationals
    @test string(rationalize(a)) == "a"
    @test string(rationalize(2a)) == "2a"
    @test string(rationalize(1.0a)) == "a"
    @test string(rationalize(0.5a)) == "(1//2)*a"
    @test string(rationalize((3 // 1) * a)) == "(3//1)*a"
    @test string(rationalize(sin(0.5a))) == "sin((1//2)*a)"
    @test string(rationalize.(v)) == "Symbolics.Num[(1//2)*a, (87//250)*b]"
    @test string(rationalize!(v)) == "Symbolics.Num[(1//2)*a, (87//250)*b]"
    @test string(v) == "Symbolics.Num[(1//2)*a, (87//250)*b]"

    # TODO rationalize doesn't rewrite a/2 since Symbolics keeps it as a division
    # with no AbstractFloat leaf for the rule to match (it used to auto-simplify
    # to 0.5a). Needs a decision: extend rationalize's match rule, or accept
    # "a / 2" as the correct current output.
    @test_broken string(rationalize.(w)) == "Symbolics.Num[(1//2)*a (87//250)*b; 3b (1//2)*a]"
    @test_broken string(rationalize!(w)) == "Symbolics.Num[(1//2)*a (87//250)*b; 3b (1//2)*a]"
    @test_broken string(w) == "Symbolics.Num[(1//2)*a (87//250)*b; 3b (1//2)*a]"

    # # Make integers
    @test string(integerize(0.0a)) == "0"
    @test string(integerize(2.0a)) == "2a"
    @test string(integerize(1.0a + a)) == "2a"
    @test string(integerize(12 // 6 * a)) == "2a"
    @test string(integerize(2.0a + 6 // 3)) == "2 + 2a"
    @test string(u) == "Symbolics.Num[2.0 + 3.0a, b]"
    @test string(integerize.(u)) == "Symbolics.Num[2 + 3a, b]"
    @test string(integerize!(u)) == "Symbolics.Num[2 + 3a, b]"
    @test string(u) == "Symbolics.Num[2 + 3a, b]"

end
