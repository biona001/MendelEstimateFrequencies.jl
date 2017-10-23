using MendelBase
using Search
using SearchSetup
using MendelEstimateFrequencies
using DataFrames 
using Distributions 

#
#This function relies on the 2 elston-stewart functions in MendelBase
#to be working properly. They do not yet have any unit testings. 
#

@testset "penetrance_estimate_frequencies" begin

end

@testset "transmission_estimate_frequencies" begin

end

@testset "prior_estimate_frequencies" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "estimate frequencies 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)

    # !locus.xlinked is true, !person.male is false
    multi_genotype = [1 0; 2 0]
    par = [0.3, 0.7]
    start = 1
    finish = 1
    i = 1
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test typeof(result) == Float64
    @test result == 0.3 * 0.7

    # !locus.xlinked is true, !person.male is true
    i = 2
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test result == 0.3 * 0.7

    # !locus.xlinked is false, !person.male is true
    locus.xlinked[:] = true
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test result == 0.3 * 0.7

    # !locus.xlinked is false, !person.male is false
    i = 3
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test result == 0.3

    # another test
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    process_keywords!(keyword, "estimate frequencies 1 Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)

    multi_genotype = [2 0; 1 0]
    par = [0.2, 0.8]
    start = 1
    finish = 1
    i = 1

    # !locus.xlinked true, !person.male false
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test typeof(result) == Float64
    @test result == 0.8 * 0.2

    # !locus.xlinked true, !person.male true
    i = 2
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test result == 0.8 * 0.2

    # !locus.xlinked false, !person.male true
    locus.xlinked[:] = true
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test result == 0.8 * 0.2

    # !locus.xlinked false, !person.male false
    i = 3
    result = MendelEstimateFrequencies.prior_estimate_frequencies(
        person, locus, multi_genotype, par, keyword, start, finish, i)
    @test result == 0.8

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

