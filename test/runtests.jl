using MendelEstimateFrequencies
using MendelBase
using Search
using SearchSetup
#
# Required external modules.
#
using DataFrames 
using Distributions 

using Base.Test

# write your own tests here

include("MendelEstimateFrequencies_test.jl")

# using Coverage
# julia -e 'Pkg.test("MendelEstimateFrequencies",coverage=true)'
# @show get_summary(process_file("src/MendelEstimateFrequencies.jl"))
