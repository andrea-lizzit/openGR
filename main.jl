using OrdinaryDiffEq
using ProgressBars
using StaticArrays
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

function f(u, p, t)
	l, θ, ϕ, pl, pt = u
	b, B2, ρ, M, a = p
	r = dneg_wh(l, ρ, M, a)
	drdl = dneg_drdl(l, ρ, M, a)
	dl = pl
	dθ = pt / r^2
	dϕ = b / (r^2 * sin(θ)^2)
	dpl = B2 * drdl / r^3
	dpt = b^2 / r^2 * cos(θ) / sin(θ)^3
	@SVector [dl, dθ, dϕ, dpl, dpt]
end


function canonical_momenta(nl, nt, np, r, θ)
	# return pl, pθ, pϕ
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

struct Bounds
	minϕ::Real
	maxϕ::Real
	minθ::Real
	maxθ::Real
end

struct WormholeParams
	ρ::Real
	M::Real
	a::Real
end

ρ = 1.
a = 0.005ρ
W = 0.10ρ
M = W / 1.42953
whparams = WormholeParams(ρ, M, a)

lc = 4.05ρ + a
θc = π/2
ϕc = 0.0;
r = dneg_wh(lc, ρ, M, a)

width, height = 500, 400; #3841, 2161; #9001, 6001;#1200, 900;


function initialize_problem(θcs, ϕcs, whparams)
	nl, np, nt = toglobal(θcs, ϕcs)
	pl, pθ, pϕ = canonical_momenta(nl, nt, np, r, θc);
	b, B2 = constants_of_motion(nl, nt, np, r, θc);
	u0 = @SVector [lc, θc, ϕc, pl, pθ];
	p = @SVector [b, B2, whparams.ρ, whparams.M, whparams.a];
	u0, p
end

function prob_func_closure(width, height, bounds, whparams)
	function prob_func_inner(prob, i, repeat)
		# problem function
		# Returns a new problem with parameters initialized from index i
		j = i % width;
		i = i ÷ width;
		θcs = bounds.minθ + (bounds.maxθ - bounds.minθ) * i / (height-1);
		ϕcs = bounds.minϕ + (bounds.maxϕ - bounds.minϕ) * j / (width-1);
		u0, p = initialize_problem(θcs, ϕcs, whparams)
		remake(prob, u0=u0, p=p)
	end
	return prob_func_inner
end

pbar = ProgressBar(total=width*height)

function output_func(sol, i)
	update(pbar);
	s = (sign(sol[end][1]) + 1) / 2
	[sol[end][2], sol[end][3], s], false
end

bounds = Bounds(0., 2π, 0., π)
p_func = prob_func_closure(width, height, bounds, whparams)
u0, p = initialize_problem(bounds.minθ, bounds.minϕ, whparams)
prob = ODEProblem(f, u0, (0.0, -1000.0), p)
eprob = EnsembleProblem(prob; output_func=output_func, prob_func=p_func)
sim = solve(eprob, Tsit5(), EnsembleThreads(), trajectories=width*height, save_everystep=false)
# bounds = (π-0.075 - 0.075/2, π-0.075, π/2 - 0.075/2, π/2) # saturn sector with artifact
# bounds = (π-0.15, π+0.15, π/2 - 0.15, π/2+0.15) # saturn sector zoom
# bounds = (π-0.15, π+0.15, π/2 - 0.15, π/2+0.15)
# bounds = (π-0.25, π+0.25, π/2 - 0.25, π/2+0.25)


image = reshape(hcat(sim...), (3, width, height))
# write image to disk in HDF5
h5write("map.h5", "map", image)