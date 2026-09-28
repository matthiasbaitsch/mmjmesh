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

    # TODO integerize only rewrites a Rational-coefficient term when it is the
    # entire top-level expression. As soon as it is summed with anything else
    # (even a bare symbol), Postwalk/maketerm judges the recursively-cleaned
    # child "equal" to the original (isequal(6//1 * a, 6a) == true) and skips
    # rebuilding the parent Add, silently discarding the cleaned version.
    @test_broken string(integerize(Num(6 // 1) * a + b)) == "6a + b"

    # # Reduce content
    @test string(cancel(Num(3) * a)) == "3a"
    @test string(cancel(Num(12) * a / Num(8))) == "(3//2)*a"

    n1 = Num(3) * a + Num(2) * b
    @test isequal(simplify(cancel(n1 / Num(5)) - n1 / Num(5); expand=true), 0)

    n2 = Num(6) * a + Num(4) * b
    @test isequal(simplify(cancel(n2 / Num(2)) - (Num(3) * a + Num(2) * b); expand=true), 0)

    # Motivating case: repeated symbolic fraction arithmetic (as done internally by
    # `integrate`) can leave numerator and denominator with a huge common factor that
    # Symbolics' own `simplify`/`expand` do not cancel.
    huge = (
        Num(202925802399989760000000) * a^4 +
        Num(131121287704608768000000) * a^2 * b^2 +
        Num(202925802399989760000000) * b^4
    ) / (Num(45528224897433600000000) * a^3 * b^3)
    expected = (Num(780) * a^4 + Num(504) * a^2 * b^2 + Num(780) * b^4) / (Num(175) * a^3 * b^3)
    @test isequal(simplify(cancel(huge) - expected; expand=true), 0)

end
