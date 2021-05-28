#!/usr/bin/env julia

using neten

adjm = data("butterfly")

netw = enhance(adjm)
