using Oceananigans
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using SeawaterPolynomials
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Oceananigans.BuoyancyModels: buoyancy_frequency
using Printf
using GLMakie

arch = CPU()

#Nx = Ny = 256
#Nz = 128

Nx = Ny = 64
Nz = 64

Lz = 200
Lx = Ly = 400

grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       halo = (5, 5, 5),
                       x = (0, Lx),
                       y = (0, Ly),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

# Representative values
T₀ = 10 # ᵒC
S₀ = 32 # g kg⁻¹
teos10_eos = TEOS10EquationOfState(; reference_density)
α = SeawaterPolynomials.thermal_expansion(T₀, S₀, 0, teos10_eos)
β = SeawaterPolynomials.haline_contraction(T₀, S₀, 0, teos10_eos)
reference_density = 1020.0
thermal_expansion = α
haline_contraction = β
linear_eos = LinearEquationOfState(; thermal_expansion, haline_contraction)
@show linear_eos
linear_buoyancy = SeawaterBuoyancy(equation_of_state=linear_eos) 
teos10_buoyancy = SeawaterBuoyancy(equation_of_state=teos10_eos) 
coriolis = FPlane(latitude=50)
advection = WENO(order=5)

Jᵘ = Field{Nothing, Nothing, Nothing}(grid)
Jᵛ = Field{Nothing, Nothing, Nothing}(grid)
Jᵀ = Field{Nothing, Nothing, Nothing}(grid)
Jˢ = Field{Nothing, Nothing, Nothing}(grid)

u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ))
v_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ))
T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵀ))
S_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jˢ))

boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs)

model = NonhydrostaticModel(; grid,
                            coriolis,
                            advection,
                            boundary_conditions,
                            buoyancy = teos10_buoyancy,
                            #buoyancy = linear_buoyancy,
                            timestepper = :RungeKutta3,
                            tracers = (:T, :S))

profiles_filename = "ocean_station_papa_profiles.jld2"
Tdata = FieldTimeSeries(profiles_filename, "T", backend=OnDisk())
Sdata = FieldTimeSeries(profiles_filename, "S", backend=OnDisk())

nᵢ = 224
Tdataᵢ = Tdata[nᵢ]
Sdataᵢ = Sdata[nᵢ]
Niz = size(Tᵢ, 3)
ic_grid = RectilinearGrid(size=(1, 1, Niz),
                          x = (0, Lx),
                          y = (0, Ly),
                          z = (-Lz, 0),
                          topology = (Periodic, Periodic, Bounded))
Tᵢ = CenterField(ic_grid)
Sᵢ = CenterField(ic_grid)
interior(Tᵢ) .= interior(Tdataᵢ)
interior(Sᵢ) .= interior(Sdataᵢ)
interior(Sᵢ, :, :, Niz) .= interior(Sᵢ, :, :, Niz-1)
fill_halo_regions!(Tᵢ)
fill_halo_regions!(Sᵢ)

T, S = model.tracers
interpolate!(T, Tᵢ)
interpolate!(S, Sᵢ)
interior(S, :, :, Nz) .= interior(S, :, :, Nz-1)

ϵ(x, y, z) = 1e-6 * (2rand() - 1)
set!(model, u=ϵ, v=ϵ, w=ϵ)

N²linear = Field(buoyancy_frequency(linear_buoyancy, model.grid, model.tracers))
N²teos10 = Field(buoyancy_frequency(teos10_buoyancy, model.grid, model.tracers))
Tav = Field(Average(T, dims=(1, 2)))
Sav = Field(Average(S, dims=(1, 2)))

N²lav = Field(Average(N²linear, dims=(1, 2)))
N²tav = Field(Average(N²teos10, dims=(1, 2)))

function plot_buoyancy_structure()
    compute!(N²lav)
    compute!(N²tav)
    compute!(Tav)
    compute!(Sav)
    fig = Figure(resolution=(1200, 600))
    axT = Axis(fig[1, 1])
    axS = Axis(fig[1, 2])
    axN = Axis(fig[1, 3])
    Tn = interior(Tav, 1, 1, :)
    Sn = interior(Sav, 1, 1, :)
    Nln = interior(N²lav, 1, 1, :)
    Ntn = interior(N²tav, 1, 1, :)
    zc = znodes(Tav)
    zf = znodes(N²av)
    scatterlines!(axT, Tn, zc)
    scatterlines!(axS, Sn, zc)
    scatterlines!(axN, Nln, zf)
    scatterlines!(axN, Ntn, zf)
    display(fig)
end

function make_heatmaps()
    u, v, w = model.velocities
    fig = Figure(resolution=(1800, 600))
    Tn = interior(T, 1, :, :)
    Sn = interior(S, 1, :, :)
    wn = interior(w, 1, :, :)
    axT = Axis(fig[1, 1])
    axS = Axis(fig[1, 2])
    axw = Axis(fig[1, 3])
    heatmap!(axT, Tn)
    heatmap!(axS, Sn)
    heatmap!(axw, wn)
    display(fig)
end

plot_buoyancy_structure()

simulation = Simulation(model, Δt=1e-2, stop_iteration=10)

function progress(sim)
    msg1 = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    u, v, w = sim.model.velocities
    msg2 = @sprintf(", max|u|: %.2e, %.2e, %.2e",
                    maximum(abs, u),
                    maximum(abs, v),
                    maximum(abs, w))

    T, S = sim.model.tracers
    msg3 = @sprintf(", extrema(T): (%.2f, %.2f)", extrema(T)...)
    msg4 = @sprintf(", extrema(S): (%.2f, %.2f)", extrema(S)...)

    @info msg1 * msg2 * msg3 * msg4
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

run!(simulation)
make_heatmaps()
plot_buoyancy_structure()

#=
fluxes_filename = "ocean_station_papa_fluxes.jld2"
τxt = FieldTimeSeries(fluxes_filename, "τx")
τyt = FieldTimeSeries(fluxes_filename, "τy")
Qt = FieldTimeSeries(fluxes_filename, "Q")
Ft = FieldTimeSeries(fluxes_filename, "F")

function set_fluxes!(sim)
    t = Time(time(sim))
    ρᵣ = reference_density
    cₚ = 3991.0

    grid = sim.model.grid
    Nz = size(grid, 3)
    S̄ = mean(interior(S, :, :, Nz))

    @inbounds begin
        τx = τxt[1, 1, 1, t]
        τy = τyt[1, 1, 1, t]
        Q  = Qt[1, 1, 1, t]
        F  = Ft[1, 1, 1, t]

        Jᵘ[1, 1, 1] = τx / ρᵣ
        Jᵛ[1, 1, 1] = τy / ρᵣ
        Jᵀ[1, 1, 1] = Q / (ρᵣ * cₚ)
        Jˢ[1, 1, 1] = S̄ * F
    end

    return nothing
end
=#

