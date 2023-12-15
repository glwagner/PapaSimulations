# This script saves temperature, salinity, and surface flux data to
# files in the Oceananigans FieldTimeSeries format, so that it can be used to
# initialize and force large eddy and single column simulations.

using ClimaOcean
using Oceananigans
using Oceananigans.Fields: interpolate!
using NCDatasets
using GLMakie

# Paths to NetCDF data
temperature_filename    = "ospapa_temperature.nc"
salinity_filename       = "ospapa_salinity.nc"
surface_fluxes_filename = "ospapa_surface_forcing.nc"
stokes_drift_filename   = "ospapa_stokes_drift.nc"

# Filenames for saved data in Oceananigans FieldTimeSeries format
profile_filename = "ocean_station_papa_profiles.jld2"
flux_filename    = "ocean_station_papa_fluxes.jld2"

isfile(profile_filename) && rm(profile_filename)

temperature_ds    = Dataset(temperature_filename)
surface_fluxes_ds = Dataset(surface_fluxes_filename)
stokes_drift_ds   = Dataset(stokes_drift_filename)
salinity_ds       = Dataset(salinity_filename)

T  = temperature_ds["T"][:, :]
S  = salinity_ds["S"][:, :]
zT = temperature_ds["z"][:]
zS = salinity_ds["z"][:]

T = reverse(T, dims=1)
S = reverse(S, dims=1)
zT = reverse(zT)
zS = reverse(zS)

tFd = surface_fluxes_ds["time"][:]
Nt = length(tFd)
hour = 3600.0
dtF = hour
tFf = (Nt - 1) * dtF
tF = collect(0:dtF:tFf)
NtF = length(tF)

tTd = temperature_ds["time"][:]
tSd = salinity_ds["time"][:]
NtT = length(tTd)
day = 24hour
dtT = day
tTf = (Nt - 1) * dtT
tT = collect(0:dtT:tTf)

mm_hour = 1e-3 / hour

u10 = surface_fluxes_ds["U10"][:]
v10 = surface_fluxes_ds["V10"][:]
τx = - surface_fluxes_ds["taux"][:]
τy = - surface_fluxes_ds["tauy"][:]
Qh = - surface_fluxes_ds["Qh"][:]
Qs = - surface_fluxes_ds["Qs"][:]
Q = Qh .+ Qs
F = - mm_hour .* surface_fluxes_ds["EMP"][:]

function make_grid(zc)
    Δz = zc[2:end] .- zc[1:end-1]
    zf = (zc[1:end-1] .+ zc[2:end]) ./ 2
    push!(zf, 0)
    pushfirst!(zf, zf[1] - Δz[1])
    Nz = length(zc)
    grid = RectilinearGrid(size=Nz, z=zf, topology=(Flat, Flat, Bounded))
    return grid
end

grid_T = make_grid(zT)
grid_S = make_grid(zS)

function field_from_data!(field, data, k, n)
    Nz = size(data, 1)
    if ismissing(data[k, n])
        if k == Nz
            kk = Nz - 1

            while kk > 1 && ismissing(data[kk, n])
                kk -= 1
            end

            if !ismissing(data[kk, n])
                field[1, 1, k] = data[kk, n]
            else
                field[1, 1, k] = 0
            end
        else
            field[1, 1, k] = field[1, 1, k+1]
        end
    else
        field[1, 1, k] = data[k, n]
    end

    return nothing
end

grid = RectilinearGrid(size=200, z=(-200, 0), topology=(Flat, Flat, Bounded))

c_loc  = (Center, Center, Center)
Tt = FieldTimeSeries(c_loc, grid, tT, backend=OnDisk(), path=profile_filename, name="T")
St = FieldTimeSeries(c_loc, grid, tT, backend=OnDisk(), path=profile_filename, name="S")

Tdata = CenterField(grid_T)
Sdata = CenterField(grid_S)

Tfine = CenterField(grid)
Sfine = CenterField(grid)

fig = Figure()
axT = Axis(fig[1, 1])
axS = Axis(fig[1, 2])

for n = 1:NtT
    local Nz
    Nz = size(T, 1)
    for k = Nz:-1:1
        field_from_data!(Tdata, T, k, n)
    end

    Nz = size(S, 1)
    for k = Nz:-1:1
        field_from_data!(Sdata, S, k, n)
    end

    fill_halo_regions!(Tdata)
    fill_halo_regions!(Sdata)

    interpolate!(Tfine, Tdata)
    interpolate!(Sfine, Sdata)

    fill_halo_regions!(Tfine)
    fill_halo_regions!(Sfine)

    @info "Saved profiles at time point $n of $NtT"

    #=
    zTdata = znodes(Tdata)
    zSdata = znodes(Sdata)
    zfine = znodes(Tfine)
    scatterlines!(axT, interior(Tfine, 1, 1, :), zfine)
    scatterlines!(axS, interior(Sfine, 1, 1, :), zfine)
    display(fig)
    sleep(1.1)
    =#

    set!(Tt, Tfine, n)
    set!(St, Sfine, n)
end

Ttmem = FieldTimeSeries(profile_filename, "T")
Stmem = FieldTimeSeries(profile_filename, "S")
t = Ttmem.times
Nt = length(t)
z = znodes(Ttmem)
fig = Figure()
axT = Axis(fig[1, 1])
axS = Axis(fig[1, 2])
slider = Slider(fig[2, 1:2], range=1:Nt, startvalue=1)
n = slider.value
Tn = @lift interior(Ttmem[$n], 1, 1, :)
Sn = @lift interior(Stmem[$n], 1, 1, :)
scatterlines!(axT, Tn, z)
scatterlines!(axS, Sn, z)
display(fig)

Q_loc  = (Center, Center, Nothing)
τxt  = FieldTimeSeries(Q_loc, grid, tF, backend=OnDisk(), path=flux_filename, name="τx")
τyt  = FieldTimeSeries(Q_loc, grid, tF, backend=OnDisk(), path=flux_filename, name="τy")
Qt   = FieldTimeSeries(Q_loc, grid, tF, backend=OnDisk(), path=flux_filename, name="Q")
Qswt = FieldTimeSeries(Q_loc, grid, tF, backend=OnDisk(), path=flux_filename, name="Qsw")
Qht  = FieldTimeSeries(Q_loc, grid, tF, backend=OnDisk(), path=flux_filename, name="Qh")
Ft   = FieldTimeSeries(Q_loc, grid, tF, backend=OnDisk(), path=flux_filename, name="F")

flux = Field{Center, Center, Nothing}(grid)

function set_flux!(flux, data, n)
    @inbounds begin
        if ismissing(data[n])
            if n > 1
                data[n] = data[n-1]
            else
                data[n] = data[n+1]
            end
        end
        flux[1, 1, 1] = data[n]
    end
    return nothing
end

for n = 1:NtF
    set_flux!(flux, τx, n)
    set!(τxt, flux, n)

    set_flux!(flux, τy, n)
    set!(τyt, flux, n)

    set_flux!(flux, Q, n)
    set!(Qt, flux, n)

    set_flux!(flux, Qs, n)
    set!(Qswt, flux, n)

    set_flux!(flux, Qh, n)
    set!(Qht, flux, n)

    set_flux!(flux, F, n)
    set!(Ft, flux, n)

    @info "Saved fluxes at time point $n of $NtF"
end

#=
fig = Figure(resolution=(3000, 2000))
axT = Axis(fig[1, 1], ylabel="z (m)", xlabel="Temperature (ᵒC)")
axS = Axis(fig[1, 2], ylabel="z (m)", xlabel="Salinity (g/kg)")

axτ = Axis(fig[2, 1:2], xlabel="Time (days)", ylabel="Momentum fluxes (N m⁻²)")
axQ = Axis(fig[3, 1:2], xlabel="Time (days)", ylabel="Heat fluxes (W m⁻²)")
axF = Axis(fig[4, 1:2], xlabel="Time (days)", ylabel="Freshwater forcing (m s⁻¹)")

linkxaxes!(axτ, axQ)
linkxaxes!(axτ, axF)

slider = Slider(fig[5, 1:2], range=1:NtT, startvalue=1)
n = slider.value
Tn = @lift T[:, $n]
Sn = @lift S[:, $n]
tn = @lift tT[$n]

scatterlines!(axT, Tn, zT)
scatterlines!(axS, Sn, zS)

xlims!(axT, 3, 14)
xlims!(axS, 32.4, 34)

lines!(axτ, ta, τx, label="τx")
lines!(axτ, ta, τy, label="τy")
vlines!(axτ, tn)
axislegend(axτ)

lines!(axQ, ta, Qh, label="Non-solar heat flux")
lines!(axQ, ta, Qs, label="Solar insolation")
lines!(axQ, ta, Q, label="Total", linewidth=3)
vlines!(axQ, tn)
axislegend(axQ)

lines!(axF, ta, F)
vlines!(axF, tn)

title = @lift string("Ocean station PAPA data at ", td[$n], " ($(n.val))")
Label(fig[0, 1:2], title)

display(fig)
=#

