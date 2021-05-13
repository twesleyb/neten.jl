#!/usr/bin/env julia

using neten
using Test
using BenchmarkTools

@testset "neten.jl" begin
    # Write your tests here.
	#hello()
	adjm = data("butterfly")
	netw = enhance(adjm)
end

# timed test
adjm = data("butterfly");
@benchmark enhance(adjm)
