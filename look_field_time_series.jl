using Oceananigans
using Oceananigans.Units
using GLMakie

profiles_filename = "ocean_station_papa_profiles.jld2"
fluxes_filename = "ocean_station_papa_fluxes.jld2"

τxt = FieldTimeSeries(fluxes_filename, "τx")
τyt = FieldTimeSeries(fluxes_filename, "τy")
Qt  = FieldTimeSeries(fluxes_filename, "Q")
Ft  = FieldTimeSeries(fluxes_filename, "F")

Tt = FieldTimeSeries(profiles_filename, "T", backend=OnDisk())
St = FieldTimeSeries(profiles_filename, "S", backend=OnDisk())

tp = Tt.times
tf = Qt.times

Ntp = length(tp)
Ntf = length(tf)

set_theme!(Theme(fontsize=36))
fig = Figure(resolution=(3000, 2000))
axT = Axis(fig[1, 1], ylabel="z (m)", xlabel="Temperature (ᵒC)")
axS = Axis(fig[1, 2], ylabel="z (m)", xlabel="Salinity (g/kg)")

axτ = Axis(fig[2, 1:2], xlabel="Time (days)", ylabel="Momentum fluxes (N m⁻²)")
axQ = Axis(fig[3, 1:2], xlabel="Time (days)", ylabel="Heat fluxes (W m⁻²)")
axF = Axis(fig[4, 1:2], xlabel="Time (days)", ylabel="Freshwater forcing (m s⁻¹)")

linkxaxes!(axτ, axQ)
linkxaxes!(axτ, axF)

slider = Slider(fig[5, 1:2], range=1:Ntp, startvalue=1)
n = slider.value

Tn = @lift interior(Tt[$n], 1, 1, :)
Sn = @lift interior(St[$n], 1, 1, :)
tn = @lift tp[$n]

scatterlines!(axT, Tn, zT)
scatterlines!(axS, Sn, zS)

xlims!(axT, 3, 14)
xlims!(axS, 32.4, 34)
ylims!(axT, -230, 30)
ylims!(axS, -230, 30)

lines!(axτ, tf, interior(τxt, 1, 1, 1, :), label="τˣ")
lines!(axτ, tf, interior(τyt, 1, 1, 1, :), label="τʸ")
vlines!(axτ, tn)
axislegend(axτ)

lines!(axQ, tf, interior(Q, 1, 1, 1, :))
vlines!(axQ, tn)

lines!(axF, tf, F)
vlines!(axF, tn)

title = @lift string("Ocean station PAPA data ",
                     tp[$n] / day, " days after March 20, 2012")
Label(fig[0, 1:2], title)

record(fig, "ocean_station_papa.mp4", 1:NtT, framerate=24) do nn
    @info "Plotting frame $n of $nn"
    n[] = nn
end

display(fig)

