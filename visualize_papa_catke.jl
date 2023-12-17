using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time
using Statistics
using GLMakie
using Printf
using Dates

profiles_filename = "ocean_station_papa_profiles.jld2"
fluxes_filename = "ocean_station_papa_fluxes.jld2"
simulation_filename = "papa_catke.jld2"

#=
τxt  = FieldTimeSeries(fluxes_filename, "τx")
τyt  = FieldTimeSeries(fluxes_filename, "τy")
Qht  = FieldTimeSeries(fluxes_filename, "Qh")
Qswt = FieldTimeSeries(fluxes_filename, "Qsw")
Ft   = FieldTimeSeries(fluxes_filename, "F")

Tobs = FieldTimeSeries(profiles_filename, "T", backend=OnDisk())
Sobs = FieldTimeSeries(profiles_filename, "S", backend=OnDisk())

ut = FieldTimeSeries(simulation_filename, "u", backend=OnDisk())
vt = FieldTimeSeries(simulation_filename, "v", backend=OnDisk())
Tt = FieldTimeSeries(simulation_filename, "T", backend=OnDisk())
St = FieldTimeSeries(simulation_filename, "S", backend=OnDisk())
κt = FieldTimeSeries(simulation_filename, "κᶜ", backend=OnDisk())
et = FieldTimeSeries(simulation_filename, "e", backend=OnDisk())
=#

t = ut.times
Nt = length(t)

n = Observable(1)
zc = znodes(ut)
zf = znodes(κt)

un = @lift interior(ut[$n], 1, 1, :)
vn = @lift interior(vt[$n], 1, 1, :)
Tn = @lift interior(Tt[$n], 1, 1, :)
Sn = @lift interior(St[$n], 1, 1, :)
κn = @lift interior(κt[$n], 1, 1, :)
en = @lift interior(et[$n], 1, 1, :)

Tnobs = @lift interior(Tobs[Time(t[$n])], 1, 1, :)
Snobs = @lift interior(Sobs[Time(t[$n])], 1, 1, :)
zobs = znodes(Tobs)

dT = 4
Navg = 10
Tmax = Ref(maximum(Tn[]) + dT)
Tavg = ones(Navg) .* Tmax[]

Tlimits = @lift begin
    if mod($n, 10) == 0
        Tn = interior(Tt[$n], 1, 1, :)
        Tavg[1:end-1] = Tavg[2:end] 
        Tavg[end] = maximum(Tn) + dT
        Tmax[] = mean(Tavg)
    end
    ((3, Tmax[]), nothing)
end

function moving_time_limits(td)
    tmin = max(-2, td - 10)
    tmax = min(td + 10, tmin + 20, td + 10, t[end] / days + 2)

    # Shift tmin back if necessary to preserve window width
    tmin = min(tmin, tmax - 20)
    return tmin, tmax
end

const τlims = (-1.2, 1.2)
const Qlims = (-1000, 400)

td = t[1] / days
τtlims = Ref(moving_time_limits(td))
Qtlims = Ref(moving_time_limits(td))

τtlimits = @lift begin
    if mod($n, 10) == 0
        td = t[$n] / days
        τtlims[] = moving_time_limits(td)
    end
    (τtlims[], τlims)
end

Qtlimits = @lift begin
    if mod($n, 10) == 0
        td = t[$n] / days
        Qtlims[] = moving_time_limits(td)
    end
    (Qtlims[], Qlims)
end

set_theme!(Theme(fontsize=36, linewidth=4))
fig = Figure(resolution=(3200, 2400))

axτ = Axis(fig[1, 1:5], limits=τtlimits, ylabel="Wind stress (N m⁻²)", xlabel="Time (days)")
axQ = Axis(fig[2, 1:5], limits=Qtlimits, ylabel="Heat flux (W m⁻²)", xlabel="Time (days)")

axu = Axis(fig[3, 1], ylabel="z (m)", xlabel="Velocities (m s⁻¹)")
axT = Axis(fig[3, 2], ylabel="z (m)", xlabel="Temperature (ᵒC)", limits=Tlimits)
axS = Axis(fig[3, 3], ylabel="z (m)", xlabel="Salinity (g kg⁻¹)")
axκ = Axis(fig[3, 4], ylabel="z (m)", xlabel="Tracer diffusivity (m² s⁻¹)")
axe = Axis(fig[3, 5], ylabel="z (m)", xlabel="Turbulent kinetic \n energy (m² s⁻²)")

#slider = Slider(fig[2, 1:5], range=1:Nt, startvalue=1)
#n = slider.value

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

tndays = @lift t[$n] / days

lines!(axτ, τxt.times ./ days, τxt[:])
lines!(axτ, τxt.times ./ days, τyt[:])
vlines!(axτ, tndays)

lines!(axQ, Qswt.times ./ days, Qswt[:])
lines!(axQ,  Qht.times ./ days,  Qht[:])
vlines!(axQ, tndays)

lines!(axu, un, zc, color=:black, label="u")
lines!(axu, vn, zc, color=:black, linestyle=:dash, label="v")

lines!(axT, Tn, zc, label="CATKE")
lines!(axT, Tnobs, zobs, linewidth=8, label="Observations")

lines!(axS, Sn, zc, label="CATKE")
lines!(axS, Snobs, zobs, linewidth=8, label="Observations")

Label(fig[0, 1:4], title)

lines!(axκ, κn, zf)
lines!(axe, en, zc)

xlims!(axu, -0.3, 0.3)
xlims!(axS, 32.4, 34)
xlims!(axe, -1e-4, 1e-2)
xlims!(axκ, -1e-2, 2e0)

axislegend(axu, position=:rb)
axislegend(axT, position=:rb)
axislegend(axS, position=:lb)

record(fig, "papa_catke.mp4", 1:Nt, framerate=96) do nn
    @info "Plotting $nn of $Nt..."
    n[] = nn
end

display(fig)

