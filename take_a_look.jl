using NCDatasets
using GLMakie

temperature_filename = "ospapa_temperature.nc"
surface_fluxes_filename = "ospapa_surface_forcing.nc"
stokes_drift_filename = "ospapa_stokes_drift.nc"
salinity_filename = "ospapa_salinity.nc"

temperature_ds    = Dataset(temperature_filename)
surface_fluxes_ds = Dataset(surface_fluxes_filename)
stokes_drift_ds   = Dataset(stokes_drift_filename)
salinity_ds       = Dataset(salinity_filename)

T = temperature_ds["T"][:, :]
S = salinity_ds["S"][:, :]
zT = temperature_ds["z"][:]
zS = salinity_ds["z"][:]

td = surface_fluxes_ds["time"][:]
Nt = length(td)
hour = 3600.0
day = 24hour
dt = hour
tf = (Nt - 1) * dt
ta = 0:dt:tf
ta = ta ./ day

td = temperature_ds["time"][:]
NtT = length(td)
tT = 0:1:(NtT-1)

hour = 3600.0
mm_hour = 1e-3 / hour

u10 = surface_fluxes_ds["U10"][:]
v10 = surface_fluxes_ds["V10"][:]
τˣ = - surface_fluxes_ds["taux"][:]
τʸ = - surface_fluxes_ds["tauy"][:]
Qh = - surface_fluxes_ds["Qh"][:]
Qs = - surface_fluxes_ds["Qs"][:]
Q = Qh .+ Qs
EP = mm_hour .* surface_fluxes_ds["EMP"][:]

set_theme!(Theme(fontsize=36))
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
ylims!(axT, -330, 30)
ylims!(axS, -330, 30)

lines!(axτ, ta, τˣ, label="τˣ")
lines!(axτ, ta, τʸ, label="τʸ")
vlines!(axτ, tn)
axislegend(axτ)

lines!(axQ, ta, Qh, label="Non-solar heat flux")
lines!(axQ, ta, Qs, label="Solar insolation")
# lines!(axQ, ta, Q, label="Total", linewidth=3)
vlines!(axQ, tn)
axislegend(axQ)

lines!(axF, ta, EP)
vlines!(axF, tn)

title = @lift string("Ocean station PAPA data at ", td[$n], " (", $n, ")")
Label(fig[0, 1:2], title)

record(fig, "ocean_station_papa.mp4", 1:NtT, framerate=24) do nn
    @info "Plotting frame $nn of $NtT"
    n[] = nn
end

display(fig)

