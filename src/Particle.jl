#=
    Particle class and constants
=#
const q = m = 1.

ArrayFloat      = Union{Vector{Float64},Vector{Vector{Float64}}}
ArrayFunction   = Union{Array{Function},Function}
# ArrayBool   = Union{Array{Bool},Bool}


struct forces
    # Set up the ODEs for the problem
    # Set up the problem
    MagneticField   :: Union{Function,Vector{Function}}
    EOM             :: Function
    EOM_GC          :: Function
    event           :: Union{Function,Nothing}
    # Event function should return the index of the field to choose
    function forces(MagField;event=nothing)
        new(MagField,MagneticForce,MagneticForce_GC,event)
    end

end


mutable struct particle
    # Position and time things
    x           :: Array{Float64}
    v           :: Array{Float64}
    gc          :: Bool
    t           :: Array{Float64}
    Δt          :: Float64
    # Functions
    B           :: Array{Float64}
    lvol        :: Vector{Int}
    gc_init     :: Bool
    x₀          :: Vector{Float64}

    # Constructor to build easier
    function particle(x::Array{Float64},v::Array{Float64},mode,Δt,Bfield,lvol;gc_initial=true)
        # Check the magnetic field type
        if typeof(Bfield) <: Function
            B = Bfield
        elseif typeof(Bfield) <: Vector{Function}
            B = Bfield[lvol]
        end

        x₀ = x
        if gc_initial
            # Convert to FO position
            x = x - guiding_center(x,v,B)
        end
        # Create the particle object
        new(x,v,mode,[0.],Δt,B(x),[lvol],gc_initial,x₀)
    end
end


mutable struct sim
    # Vector that holds particle
    nparts  :: Int64 
    sp      :: Vector{particle}

    function sim(nparts::Int,x₀::Array{Float64},v₀::Array{Float64},mode::Bool,Δt::Vector{Float64},Bfield,lvol)
        
        if typeof(x₀) == Vector{Float64}
            # Ensure position can be read
            x₀ = [x₀ for i in 1:nparts]
        end
        if typeof(v₀) == Vector{Float64}
            # Ensure velocity can be read
            v₀ = [v₀ for i in 1:nparts]
        end
        if typeof(lvol) == Int64
            # Ensure lvol can be read
            lvol = [lvol for i in 1:nparts]
        end

        parts = Array{particle}(undef,nparts)
        for i = 1:nparts
            parts[i] = particle(x₀[i],v₀[i],mode,Δt[i],Bfield,lvol[i])
        end
        new(nparts,parts)
    end
end


mutable struct exact_particle
    # Structure for hold the exact solutions to particles
    x           :: Array{Float64}
    v           :: Array{Float64}
    x_boundary  :: Array{Float64}
    v_boundary  :: Array{Float64}
    t           :: Vector{Float64}
    t_boundary  :: Vector{Float64}
end


mutable struct analytic_sim
    # Container for multiple analytic particles
    nparts  :: Int64
    sp      :: Array{exact_particle}

    function analytic_sim(nparts::Int64)
        parts = Array{exact_particle}(undef,nparts)
        new(nparts,parts)
    end
end



#==
    CONSTRUCTORS
==#




#=====
    Particle based fns
=====#

function guiding_center(x::Vector{Float64},v::Vector{Float64},Bfield::Function)
    B = Bfield(x)
    gc = m/q * cross(v,B)/norm(B,2)^2
    return gc
end

function MagneticForce_GC(xv::Array{Float64},t::Float64,B::Array)
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

