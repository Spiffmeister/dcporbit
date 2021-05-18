#=
    Particle class and constants
=#
const q = m = 1.

ArrayFloat  = Union{Array{Float64},Float64}
# ArrayBool   = Union{Array{Bool},Bool}



mutable struct particle
    # Position and time things
    x           :: Array{Float64}
    v           :: Array{Float64}
    gc          :: Bool
    t           :: Array{Float64}
    Δt          :: Float64
    # Functions
    dvdt        :: Function
    Bfield      :: Union{Array{Function},Function}
    B           :: Array{Float64}
    lvol        :: Vector{Int}
    # Constructor to build easier
    function particle(x,v,mode,Δt,dvdt,Bfield,lvol;gc_initial=true)
        typeof(Bfield) <: Function ? B = Bfield : B = Bfield[lvol]
        if gc_initial
            # Convert to FO position
            x = x - guiding_center(x,v,B)
        end
        new(x,v,mode,[0.],Δt,dvdt,Bfield,B(x),[lvol])
    end
end

mutable struct exact_particle
    x           :: Array{Float64}
    v           :: Array{Float64}
    x_boundary  :: Array{Float64}
    v_boundary  :: Array{Float64}
    t           :: Vector{Float64}
    t_boundary  :: Vector{Float64}
end

mutable struct sim
    # Vector that holds particle
    sp  :: Vector{particle}

    function sim(nparts::Int,x₀::Vector{Vector{Float64}},v₀::Vector{Vector{Float64}},mode::Bool,Δt::Vector{Float64},dvdt,Bfield,lvol)
        parts = Array{particle}(nparts)

        if typeof(x₀) == Array{Float64}
            # Ensure position can be read
            x₀ = [x₀ for i in 1:nparts]
        end
        if typeof(v₀) == Array{Float64}
            # Ensure velocity can be read
            v₀ = [v₀ for i in 1:nparts]
        end
        if typeof(lvol) == Int64
            # Ensure lvol can be read
            lvol = [lvol for i in 1:nparts]
        end

        if length(x) == nparts
            for i = 1:nparts
                # Check if Bfield is a function or an array
                parts[i] = particle(x₀[i],v₀[i],mode,Δt[i],dvdt,Bfield,lvol[i])
            end
        else

        end
        new(parts)
    end

end

mutable struct analytic
    sp  :: Array{exact_particle}
end



#==
    CONSTRUCTORS
==#




#=====
    Particle based fns
=====#

function guiding_center(x::Array{Float64},v::Array{Float64},Bfield::Function)
    B = Bfield(x)
    gc = m/q * cross(v,B)/norm(B,2)^2
    return gc
end

function ODEgc(xv::Array{Float64},t::Float64,B::Array)
    # x = xv[1:3]
    v = xv[4:6]
    v = dot(v,B)/norm(B,2)^2 * B
    dvdt = zeros(3)
    xv = vcat(v,dvdt)
    return xv
end

function MagneticForce(xv::Array{Float64},t::Float64,B::Array)
    # x = xv[1:3]
    v = xv[4:6]
    dvdt = q/m*cross(v,B)
    xv = vcat(v,dvdt)
    return xv
end

