using DataFrames
using CSV
using Plots
using Unitful
using UnitfulRecipes
using Measurements
using LsqFit
using Latexify
using UnitfulLatexify

# Set some plotting defaults for later
# Unfortunately it doesn't accept unitformat here
default(legend=:topleft, minorgrid=:true)

# Load in data
df = DataFrame(CSV.File("data/dataResistance.csv", header=2, skipto=3))

# Add appropriate units to dataframe columns
df[!, "f(kHz)"] = Quantity.(df[!, "f(kHz)"], u"kHz")
df[!, 2:8] = Quantity.(df[!, 2:8], u"Ω")

# Frequency and η = ratio of resistances Rac/Rdc
freq = df[:,1]
η(col) = df[:, col] ./ df[1, col]

# Plot η against freq
plt = plot(ylabel="\\eta", unitformat=latexify)
for col in 2:8
    plot!(freq, η(col), xlabel="f", label=names(df)[col], unitformat=latexroundunitlabel)
end
display(plt)

# Plot η against √freq
sqrtPlt = plot(ylabel="\\eta", unitformat=latexify)
for col in 2:8
    plot!(sqrt.(freq), η(col), xlabel="f^{\\frac{1}{2}}", label=names(df)[col], unitformat=latexroundunitlabel)
end
display(sqrtPlt)

# Straight line model function with unknown parameters p
# It would be nested into the fit function, but it is also used in FitPlot
m(x, p) = p[1] .* x .+ p[2]

"""
Fit straight line to straight section of log-log plot. Takes a df column number for input.
"""
function fit(col::Integer, stripbelow::Real=0.1)
    ydata = log.(η(col))
    # Now we want to strip out the inital bit
    i = findfirst(x -> x >= stripbelow, ydata)
    if i === nothing
        i = 1
    end
    ydata = ydata[i:end]
    xdata = log.(freq[i:end] ./ 1u"Hz")
    p0 = [0.5, 8] # Initial parameter guess
    fitResults = curve_fit(m, xdata, ydata, p0)
    return fitResults.param
end

@userplot FitPlot
@recipe function f(o::FitPlot)
    params = o.args
    x = LinRange(0, 14, 14 * 4)
    y = m.(x, params)
    i = findfirst(x -> x > -0.1, y)
    y = y[i:end]
    x = x[i:end]

    label := false
    linestyle := :dash

    x, y
end

xintercept(p) = - p[2] / p[1]

dfAnalysed = DataFrame(Material=names(df)[2:end], fc=0.0u"Hz")

# Plot ln(η) against ln(freq)
logPlt = plot(ylabel="ln(η)")
for col in 2:8
    plot!(log.(freq ./ 1u"Hz"), log.(η(col)), xlabel="ln(f/Hz)", label=names(df)[col], unitformat=latexroundunitlabel) # Plot data
    color = logPlt[1][end][:linecolor] # Grab colour of latest line
    params = fit(col) # Lst sq fit to straight section
    xint  = xintercept(params) # Find x-intercept
    display(xint)
    fc = exp(xint) * u"Hz" # Find critical frequency from x-intercept
    dfAnalysed[col - 1, 2] = fc # Store critical frequency
    # fcPretty = round(xintercept(params); digits=2) # Prettify it
    # logPlt[1][end][:label] = "$(logPlt[1][end][:label]) fc = $fcPretty"
    fitplot!(params, linecolor=color) # Plot fitted line in same colour as source data
end
display(logPlt)

# Load in diameters
diametersRaw = DataFrame(CSV.File("data/dataDiameters.csv"))
avgDiam = Array(diametersRaw[end - 1, 2:end - 2]) * 1u"mm"
errDiam = Array(diametersRaw[end, 2:end - 2]) * 1u"mm"
diams = avgDiam .± errDiam
# Account for μ metal sample being rectangular, not cylindrical
μWidth = (diametersRaw[end - 1,end - 1] * 1u"mm") ± (diametersRaw[end,end - 1] * 1u"mm")
μHeight = (diametersRaw[end - 1,end] * 1u"mm") ± (diametersRaw[end,end] * 1u"mm")
μDiameter = 2 * μWidth * μHeight / (μWidth + μHeight)
push!(diams, μDiameter)
# Add to dataframe
dfAnalysed[!, :Radius] = diams ./ 2

# Conductivity, from online 
σ_copper = 58.7e6u"S/m"
σ_steel = 10.1e6u"S/m"
σ = [fill(σ_copper, 4); fill(σ_steel, 2); 1/(60u"μΩ*cm")]
dfAnalysed[!, :σ] = σ

# Find μ
μ(fc, r, σ) = uconvert(u"H/m", 4 / (π * r^2 * σ * fc))
μRelative = μ.(dfAnalysed[:, :fc][1:end], dfAnalysed[:, :Radius][1:end], dfAnalysed[:, :σ][1:end]) / u"μ0"
dfAnalysed[!, :μRel] = μRelative

# Add η @ f=1MHz as column
dfAnalysed[!, :η1MHz] = Array(df[13, 2:end]) ./ Array(df[1, 2:end])
copper = sort(dfAnalysed[1:4, :], :Radius)
steel = sort(dfAnalysed[5:6, :], :Radius)


radiusPlt = plot(copper[!, :Radius], copper[!, :η1MHz], label="Copper", xlabel="Radius", ylabel="η")
plot!(steel[!, :Radius], steel[!, :η1MHz], label="Steel")
display(radiusPlt)