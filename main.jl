using DifferentialEquations
using ProgressBars
using HDF5

function dneg_wh(l, ρ, M, a)
	if (abs(l) > a)
		x = 2(abs(l) - a) / (π*M)
		return ρ + M*(x*atan(x) - 0.5log(1 + x^2))
	else
		return ρ
	end
end

function dneg_drdl(l, ρ, M, a)
	if (abs(l) > a)
		return 2/π * atan(2(abs(l)-a)/(π*M))
	else
		return 0
	end
end

function f!(du, u, p, t)
	l, θ, ϕ, pl, pt = u
	b, B2, fr, fdrdl = p
	r = fr(l)
	drdl = fdrdl(l)
	du[1] = pl
	du[2] = pt / r^2
	du[3] = b / (r^2 * sin(θ)^2)
	du[4] = B2 * drdl / r^3
	du[5] = b^2 / r^2 * cos(θ) / sin(θ)^3
end

function dummyf!(du, u, p, t)
	l, θ, ϕ, pl, pt = u
	b, B2, fr, fdrdl = p
	r = fr(l)
	drdl = fdrdl(l)
	du[1] = pl
	du[2] = pt / r^2
	du[3] = b / (r^2 * sin(θ)^2)
	du[4] = B2 * drdl / r^3
	du[5] = b^2 / r^2 * cos(θ) / sin(θ)^3
end

function canonical_momenta(nl, nt, np, r, θ)
	# return pl. pθ. pϕ
	return nl, r*nt, r*sin(θ)*np
end

function constants_of_motion(nl, nt, np, r, θ)
	# return b, B2
	return r * sin(θ) * np, r^2 * (nt^2 + np^2)
end

function toglobal(θ, ϕ)
	# transform polar coordinate direction in the local sky to global spherical polar coordinates
	return -sin(θ)*cos(ϕ), -sin(θ)*sin(ϕ), cos(θ)
end

ρ = 1.
a = 0.005ρ
W = 0.10ρ
M = W / 1.42953
lc = 4.05ρ + a
θc = π/2
ϕc = 0.0;
r = dneg_wh(lc, ρ, M, a)

width, height = 600, 400; #3841, 2161; #9001, 6001;#1200, 900;
image = zeros(height, width, 3);
bounds = (0., 2π, 0., π)
# bounds = (π-0.075 - 0.075/2, π-0.075, π/2 - 0.075/2, π/2) # saturn sector with artifact
# bounds = (π-0.15, π+0.15, π/2 - 0.15, π/2+0.15) # saturn sector zoom
# bounds = (π-0.15, π+0.15, π/2 - 0.15, π/2+0.15)
# bounds = (π-0.25, π+0.25, π/2 - 0.25, π/2+0.25)
fr = l -> dneg_wh(l, ρ, M, a)
fdrdl = l -> dneg_drdl(l, ρ, M, a)
for (i, θcs) = enumerate(ProgressBar(LinRange(bounds[3], bounds[4], height))), (j, ϕcs) = enumerate(LinRange(bounds[1], bounds[2], width))
	nl, np, nt = toglobal(θcs, ϕcs)
	pl, pθ, pϕ = canonical_momenta(nl, nt, np, r, θc)
	b, B2 = constants_of_motion(nl, nt, np, r, θc)
	u0 = [lc, θc, ϕc, pl, pθ]
	tspan = (0.0, -1000.0)
	prob = ODEProblem(f!, u0, tspan, (b, B2, fr, fdrdl))
	sol = solve(prob, reltol=1e-6)
	l, θ, ϕ = sol(-1000.0, idxs=1:3)
	# l, θ, ϕ = 1, θcs, ϕcs # use this for sanity check
	s = (sign(l) + 1) / 2
	image[i, j, :] = [θ, ϕ, s]; # [mod(θ, π), mod(ϕ, 2π), s]
end

# write image to disk in HDF5
h5write("map.h5", "map", permutedims(image, [3, 2, 1]))


# using Plots;
# heatmap(image)
# make a histogram of image[:, :, 1]
# thetas = reshape(image[:, :, 1], :)
# phis = reshape(image[:, :, 2], :)
# histogram(phis, bins=100)