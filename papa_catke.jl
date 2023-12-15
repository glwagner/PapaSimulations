using Oceananigans
using Oceananigans.Units

using Oceananigans.Grids: znode, zspacing
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.BuoyancyModels: buoyancy_frequency

using SeawaterPolynomials
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf
using GLMakie

import Dates

Δz                           = 10   # grid spacing
Lz                           = 200  # domain depth
reference_density            = 1020.0
teos10_heat_capacity         = 3991.86795711963
shortwave_penetration_length = 16   # meters
ocean_station_papa_latitude  = 50.1 # deg
seawater_albedo              = 0.08
Nz = ceil(Int, Lz / Δz) # equal spacing

grid = RectilinearGrid(size=Nz, z=(-Lz, 0), topology=(Flat, Flat, Bounded))

# Forcing and surface flux time-series
Jᵘ = Ref(0.0)
Jᵛ = Ref(0.0)
Jᵀ = Ref(0.0)
Jˢ = Ref(0.0)
I₀ = Ref(0.0)

#=
fluxes_filename = "ocean_station_papa_fluxes.jld2"
τxt  = FieldTimeSeries(fluxes_filename, "τx")
τyt  = FieldTimeSeries(fluxes_filename, "τy")
Qht  = FieldTimeSeries(fluxes_filename, "Qh")
Qswt = FieldTimeSeries(fluxes_filename, "Qsw")
Ft   = FieldTimeSeries(fluxes_filename, "F")
=#

function set_fluxes!(sim)
    t = Time(time(sim))
    ρᵣ = reference_density
    cₚ = teos10_heat_capacity
    α = seawater_albedo

    grid = sim.model.grid
    Nz = size(grid, 3)

    @inbounds begin
        τx  = τxt[1, 1, 1, t]
        τy  = τyt[1, 1, 1, t]
        Qh  = Qht[1, 1, 1, t]
        Qsw = Qswt[1, 1, 1, t]
        F   = Ft[1, 1, 1, t]
        S̃   = S[1, 1, Nz]
    end

    Jᵘ[] = τx / ρᵣ
    Jᵛ[] = τy / ρᵣ
    Jᵀ[] = Qh / (ρᵣ * cₚ)
    Jˢ[] = S̃ * F
    I₀[] = - (1 - α) * Qsw

    return nothing
end

u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵘ))
v_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵛ))
T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jᵀ))
S_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jˢ))

# I = I₀ exp(z / λ)
# → δₖI / Δz = I₀ (exp(z⁺ / λ) - exp(z⁻ / λ)) / Δz
@inline function shortwave_flux_divergence(i, j, k, grid, clock, fields, p)
    z⁺ = znode(i, j, k+1, grid, Center(), Center(), Face())
    z⁻ = znode(i, j, k,   grid, Center(), Center(), Face())
    δI = p.I₀[] * (exp(z⁺ / p.λ) - exp(z⁻ / p.λ))

    Δz = zspacing(i, j, k, grid, Center(), Center(), Center())
    return δI / Δz / (p.ρᵣ * p.cₚ)
end

parameters = (I₀ = I₀,
              λ  = shortwave_penetration_length,
              ρᵣ = reference_density,
              cₚ = teos10_heat_capacity)

solar_insolation = Forcing(shortwave_flux_divergence; parameters, discrete_form=true)

equation_of_state = TEOS10EquationOfState(; reference_density)
buoyancy = SeawaterBuoyancy(; equation_of_state)
coriolis = FPlane(latitude=ocean_station_papa_latitude)

boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs)
closure = CATKEVerticalDiffusivity()

model = HydrostaticFreeSurfaceModel(; grid,
                                    coriolis,
                                    boundary_conditions,
                                    closure,
                                    buoyancy,
                                    forcing = (; T=solar_insolation),
                                    tracers = (:T, :S, :e))

profiles_filename = "ocean_station_papa_profiles.jld2"
Tdata = FieldTimeSeries(profiles_filename, "T", backend=OnDisk())
Sdata = FieldTimeSeries(profiles_filename, "S", backend=OnDisk())

nᵢ = 100
T, S, e = model.tracers
interpolate!(T, Tdata[nᵢ])
interpolate!(S, Sdata[nᵢ])
set!(model, e=1e-6)
model.clock.time = Tdata.times[nᵢ]
stop_time = 365days

simulation = Simulation(model; Δt=10minutes, stop_time)

function progress(sim)
    msg = string("Iter: ", iteration(sim), ", time: ", prettytime(sim))

    T, S, e = sim.model.tracers
    κᶜ = model.diffusivity_fields.κᶜ
    msg *= @sprintf(", extrema(T): (%.2f, %.2f)", extrema(T)...)
    msg *= @sprintf(", extrema(S): (%.2f, %.2f)", extrema(S)...)
    msg *= @sprintf(", max e: %.2e", maximum(e))
    msg *= @sprintf(", max κ: %.2e", maximum(κᶜ))

    @info msg
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
simulation.callbacks[:fluxes] = Callback(set_fluxes!)

κᶜ = model.diffusivity_fields.κᶜ
κᵘ = model.diffusivity_fields.κᵘ
outputs = merge(model.velocities, model.tracers, (; κᶜ, κᵘ))
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(30minutes),
                                                    filename = "papa_catke.jld2",
                                                    overwrite_existing = true)

run!(simulation)

Tobs = FieldTimeSeries(profiles_filename, "T", backend=OnDisk())
Sobs = FieldTimeSeries(profiles_filename, "S", backend=OnDisk())

filename = "papa_catke.jld2"
ut = FieldTimeSeries(filename, "u", backend=OnDisk())
vt = FieldTimeSeries(filename, "v", backend=OnDisk())
Tt = FieldTimeSeries(filename, "T", backend=OnDisk())
St = FieldTimeSeries(filename, "S", backend=OnDisk())
κt = FieldTimeSeries(filename, "κᶜ", backend=OnDisk())
et = FieldTimeSeries(filename, "e", backend=OnDisk())
t = ut.times
Nt = length(t)

set_theme!(Theme(fontsize=36))
fig = Figure(resolution=(2000, 800))
axu = Axis(fig[1, 1], ylabel="z (m)", xlabel="Velocities (m s⁻¹)")
axT = Axis(fig[1, 2], ylabel="z (m)", xlabel="Temperature (ᵒC)")
axS = Axis(fig[1, 3], ylabel="z (m)", xlabel="Salinity (g kg⁻¹)")
axκ = Axis(fig[1, 4], ylabel="z (m)", xlabel="Tracer diffusivity (m² s⁻¹)")
axe = Axis(fig[1, 5], ylabel="z (m)", xlabel="Turbulent kinetic \n energy (m² s⁻²)")
slider = Slider(fig[2, 1:4], range=1:Nt, startvalue=1)
n = slider.value
z = znodes(ut)

un = @lift interior(ut[$n], 1, 1, :)
vn = @lift interior(vt[$n], 1, 1, :)
Tn = @lift interior(Tt[$n], 1, 1, :)
Sn = @lift interior(St[$n], 1, 1, :)
κn = @lift interior(κt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)

Tnobs = @lift interior(Tobs[Time(t[$n])], 1, 1, :)
Snobs = @lift interior(Sobs[Time(t[$n])], 1, 1, :)
zobs = znodes(Tobs)

lines!(axu, un, z, color=:black, label="u")
lines!(axu, vn, z, color=:black, linestyle=:dash, label="v")

lines!(axT, Tn, z, label="CATKE")
lines!(axT, Tnobs, zobs, linewidth=5, label="Observations")

lines!(axS, Sn, z, label="CATKE")
lines!(axS, Snobs, zobs, linewidth=5, label="Observations")

zκ = znodes(κt)
lines!(axκ, κn, zκ)
lines!(axe, en, z)

title = @lift begin
    date = Dates.DateTime(2012, 3, 20) + Dates.Second(round(Int, t[$n]))
    yr = Dates.year(date)
    mo = Dates.monthname(date)
    dy = Dates.day(date)
    hr = Dates.hour(date)
    min = Dates.minute(date)

    string("Ocean Station Papa observations vs CATKE predictions",
           " on $mo $dy, $yr at ", @sprintf("%02d:%02d", hr, min))
end
           
Label(fig[0, 1:4], title)

xlims!(axu, -0.5, 0.5)
xlims!(axT, 3, 12)
xlims!(axS, 32.4, 34)
xlims!(axe, -1e-4, 1e-3)
xlims!(axκ, -1e-4, 1e-1)
axislegend(axu, position=:rb)
axislegend(axT, position=:rb)
axislegend(axS, position=:lb)

display(fig)
