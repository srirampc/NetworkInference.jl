if size(ARGS, 1) < 2
    println("Not enough Arguments :: ", ARGS)
    println("Usage :: runPIDC.jl /path/to/input_file /path/to/output_file ngenes [Pkg Environment Location] ")
    exit(1)
end

println("Running with Input Args :: ", ARGS)
# if given, third argument is environment
if size(ARGS, 1) > 3
    env_name = string(ARGS[4])
    import Pkg
    println("Activating Environment :: ", env_name)
    Pkg.activate(env_name)
end
# import packages
import NetworkInference
import LightGraphs

dataset_file = string(ARGS[1])
output_file = string(ARGS[2])
num_genes = parse(Int32, string(ARGS[3]))

# load input file
if endswith(dataset_file, ".h5ad")
    println("Loading h5ad dataset :: ", dataset_file,
            " w. ", num_genes, " genes." )
    @time genes = NetworkInference.get_h5ad_nodes(dataset_file,
                                                  number_of_nodes=num_genes);
else
    println("Loading dataset :: ", dataset_file)
    @time genes = NetworkInference.get_nodes(dataset_file);
end
 
# generate network
println("Building Network h5ad dataset :: ", dataset_file)
algorithm = NetworkInference.PIDCNetworkInference()
@time network = NetworkInference.InferredNetwork(algorithm, genes);
 

# write output file 
println("Writing Network to file :: ", output_file)
NetworkInference.write_network_file(output_file, network);
