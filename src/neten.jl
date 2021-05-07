#!/usr/bin/env julia
#
# title: neten.jl
# author: twab
# description: julia implementation of network enhancement
#
module neten

## ---- imports

using RData
using LinearAlgebra



## ---- functions

function load_rdata(rdata::String, package::String)
	# FIXME: why not have function data(String)?
	# load any data object associated with an R package
	# NOTE: this doesn't work as expected bc most package data is not saved
	# internally as .rda
	# construct path to R data.rda file in data/
	file = joinpath("data","$rdata.rda")
	# construct R command
	cmd = "system.file('$file', package='$package')"
	# parse output of R command
	out = read(`Rscript -e $cmd`,String)
	path = strip(split(out)[2],'\"')
	return load(path)[rdata]
end


function colSums(a::Array{Float64,2}, dropdim::Bool=true)
	# return column sums of an array
	col_sums = sum(a, dims=1)
	if dropdim
		col_sums = dropdims(col_sums, dims=1)
	end
	return col_sums
end


function rowSums(a::Array{Float64,2}, dropdim::Bool=true)
	# return row sums of an array
	row_sums = sum(a, dims=2)
	if dropdim
		row_sums = dropdims(row_sums, dims=2)
	end
	return row_sums
end


function NE_dn(Win::Array{Float64,2}, k::Int64=size(Win,1)) 
	# what is this function doing? 
	# NE = network enhancment?, dn == ??
	#k = size(Win, 1)
	Wk = Win * k
	d = rowSums(abs.(Wk))
	D = diagm(1 ./ d)
	Wne = D * Wk
	return Wne
end


function dominateset(Win::Array{Float64,2}, k::Int64=minimum([20,size(Win,1)-1]))
	# dominating set is set of connected nodes
	W0 = abs.(Win)
	#k = minimum([k, size(W0, 1)-1])
	# sort each row in decreasing order
	# do this to each row!
	A = sort(W0, dims=2, rev=true)
	B = mapslices(x -> sortperm(x, rev=true), 
		      W0; dims=2)
	res = A[:,1:k]
	x = collect(1:size(W0, 1))
	inds=repeat(x, outer = [1, k])
	loc = B[:, 1:k]
	P0 = zeros(size(W0,2), size(W0,1))
	idx = (vec(loc .- 1)) * size(W0,2) + vec(inds)
	P0[idx] = res
	P = (P0 + transpose(P0)) / 2
	return P .* sign.(Win)
end


function transitionFields(Win::Array{Float64,2})
	zeroindex = rowSums(Win).==0
	W1 = Win * size(Win,1)
	W2 = NE_dn(W1)
	# do this operation on each col / sqrt(sum(abs(col))
	W3 = mapslices(x -> x./sqrt(sum(abs.(x)) + eps()), 
		       W2; dims=1)
	# compute cross product: x %*% t(x) (tcrossprod)
	W4 = W3 * transpose(W3)
	W4[zeroindex,:] .= 0
	W4[:,zeroindex] .= 0
	return W4
end


function neten(Win::Array{Float64,2}, k::Int64=20, alpha::Float64=0.9, diffusion::Float64=2.0)
	# network enhancement diffusion process is given by:
	# W_{t+1} = αΓ x W_{t} x Γ + Γ(1 - α) (eq 2)
	# set diag == 0
	Win[diagind(Win)] .= 0;
	# remove unconnected nodes  
	nonzero = colSums(Win) .> 0;
	W0 = Win[nonzero, nonzero];
	# do NE_dn
	Wne = NE_dn(W0);
	# mean? or some sort of normalization
	Wnorm = (Wne + transpose(Wne)) ./ 2;
	# calc dominating set
	P = dominateset(Wnorm, k);
	# working with P...
	I = diagm(ones(size(P,1))); # identiy matrix
	X = diagm(rowSums(abs.(P))); # diagonal = rowSums
	Pnew = P .+ I .+ X;
	# do transitionFields 
	Pt = transitionFields(Pnew)
	# compute eigen vectors and values
	e = eigen(Pt)
	d = e.values[end:-1:1] .- eps() # rev order to match Rout
	U = -1 .* e.vectors[:, end:-1:1] # rev order and flip sign to match R output
	# generate diag matrix D from difusion process
	dd = (1-alpha) * d ./ (1 .- alpha * d.^diffusion) # <- diffusion process
	D = diagm(dd)
	W2 = U * D * transpose(U) 
	W3 = W2 .* (1 .- diagm(ones(size(W2,1)))) ./ (1 .- diag(W2))
	# DD is is the diagonal matrix whose entries are the degree of the vertices
	DD = diagm(colSums(abs.(W0)))
	W4 = DD * W3 
	# remove negative
	W4[W4 .< 0] .= 0
	W5 = (W4 + transpose(W4)) / 2
	Wout = zeros(size(Win))
	Wout[nonzero, nonzero]  = W5
	# set diag to zero
	#Wout[diagind(Wout)] .= 0
	return Wout
end

## ---- exports

@export neten

end
