#=
    Particle class and constants
=#
const q = m = 1.

mutable struct particle
    # Position and time things
    x   :: Array{Float64}
    v   :: Array{Float64}
    gc  :: Bool
    t   :: Array{Float64}
    Δt  :: Float64
    # Functions
    dvdt :: Function
    Bfield :: Union{Array{Function},Function}
    B :: Array{Float64}
    lvol :: Vector{Int}
    # Constructor to build easier
    function particle(x,v,mode,t,Δt,dvdt,Bfield,lvol;gc_initial=true)
        typeof(Bfield) <: Function ? B = Bfield : B = Bfield[lvol]
        if gc_initial
            # Convert to FO position
            x = x - guiding_center(x,v,B)
        end
        new(x,v,mode,[t],Δt,dvdt,Bfield,B(x),[lvol])
    end
end


#=
    Particle based fns
=#

function guiding_center(x::Array{Float64},v::Array{Float64},Bfield::Function)
    B = Bfield(x)
    gc = m/q * cross(v,B)/norm(B,2)^2
end

function ODEgc(xv::Array{Float64},t::Float64,B::Array)
    x = xv[1:3]
    v = xv[4:6]
    v = dot(v,B)/norm(B,2)^2 * B
    dvdt = zeros(3)
    xv = vcat(v,dvdt)
    return xv
end

function MagneticForce(xv::Array{Float64},t::Float64,B::Array)
    x = xv[1:3]
    v = xv[4:6]
    dvdt = q/m*cross(v,B)
    xv = vcat(v,dvdt)
    return xv
end

