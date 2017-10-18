using MendelBase
using Search
using SearchSetup
using MendelEstimateFrequencies
#
# Required external modules.
#
using DataFrames 
using Distributions 

@testset "basics" begin
    # keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    # process_keywords!(keyword, "estimate frequencies 1 Control.txt", "")
    # (pedigree, person, nuclear_family, locus, snpdata,
    # locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    # read_external_data_files(keyword)
    # skipped_loci = 0

    # keyword["eliminate_genotypes"] = true

    # locus_frame[:PedFrequency] = zeros(size(locus_frame, 1))
    # n = 0

    # (model_loci, locus.model_loci) = (locus.model_loci, 1)
    # model_locus = similar(locus.model_locus)
    # copy!(model_locus, locus.model_locus)
    # locus.model_locus = zeros(Int, 1)


    # loc = 1 
    # locus.model_locus[1] = loc
    # keyword["constraints"] = 1
    # keyword["goal"] = "maximize"
    # keyword["parameters"] = locus.alleles[loc]
    # keyword["title"] = "Estimate Frequencies analyis for " * locus.name[loc]
    # parameter = set_parameter_defaults(keyword)
    # parameter =
    #      MendelEstimateFrequencies.
    #      initialize_optimization_estimate_frequencies!(locus, parameter, keyword)
end

@testset "penetrance_estimate_frequencies" begin
    # function not written yet?
end

@testset "transmission_estimate_frequencies" begin
    
end

@testset "initialize_optimization_estimate_frequencies" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "estimate frequencies 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)

    # keyword["eliminate_genotypes"] = true
    # locus_frame[:PedFrequency] = zeros(size(locus_frame, 1))

    (model_loci, locus.model_loci) = (locus.model_loci, 1)
    model_locus = similar(locus.model_locus)
    copy!(model_locus, locus.model_locus)
    locus.model_locus = zeros(Int, 1)
    locus.model_locus[1] = 1
    keyword["constraints"] = 1
    keyword["goal"] = "maximize"
    keyword["parameters"] = locus.alleles[1]
    keyword["title"] = "Estimate Frequencies analyis for " * locus.name[1]
    parameter = set_parameter_defaults(keyword)
    parameter = MendelEstimateFrequencies.
         initialize_optimization_estimate_frequencies!(
            locus, parameter, keyword)

    @test parameter.cases == 0
    @test parameter.constraints == 1
    @test parameter.goal == "maximize"
    @test typeof(parameter.output_unit) <: IO 
    @test parameter.output_unit.name == 
        "<file estimate frequencies 1 Output.txt>"
    @test parameter.parameters == 2
    @test parameter.points == 0
    @test parameter.standard_errors == false
    @test parameter.title == "Estimate Frequencies analyis for SNP1"
    @test parameter.travel == "search"
    @test size(parameter.constraint) == (1, 2)
    @test eltype(parameter.constraint) == Float64
    @test parameter.constraint[1, 1] == 1.0
    @test parameter.constraint[1, 2] == 1.0
    @test parameter.constraint_level[1] == 1.0
    @test parameter.function_value[1] == 0.0
    @test parameter.function_value[2] == 0.0
    @test size(parameter.grid) == (0, 2) #not using grid
    @test parameter.min[1] == 1.0e-5
    @test parameter.min[2] == 1.0e-5
    @test parameter.max[1] == Inf
    @test parameter.max[2] == Inf
    @test typeof(parameter.name) <: Array{AbstractString}
    @test parameter.par[1] == 0.6811
    @test parameter.par[2] == 0.3189
end

