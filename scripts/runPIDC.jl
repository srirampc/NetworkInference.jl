using ArgParse

function parse_commandline()
    argpv = ArgParseSettings()

    @add_arg_table argpv begin
        "--env", "-e"
            help = "Path to Environment"
        "--hdf5"
            help = "Use HDF5 interface"
            action = :store_true
        "--transpose"
            help = "Transpose input matrix"
            action = :store_true
        "--num_genes", "-n"
            help = "Number of Genes"
            arg_type = Int
            default = 0
        "--genes_path", "-g"
            help = "HDF5 dataset path to list of gene identifiers"
            default = "/var/gene_ids"
        "--data_path", "-d"
            help = "HDF5 dataset path to datasets"
            default = "/X"
        "input_file"
            help = "Path to input anndata/HDF5/csv file"
            required = true
        "output_file"
            help = "Path to output csv file"
            required = true
    end

    return parse_args(argpv)
end

input_args = parse_commandline()
println("Running with Input Args :: ", input_args)

dataset_file = input_args["input_file"]
output_file = input_args["output_file"]
data_path = input_args["data_path"]
var_path = input_args["genes_path"]
num_genes =  input_args["num_genes"]
h5_flag = input_args["hdf5"]

# if given, third argument is environment
if !isnothing(input_args["env"])
    env_name = input_args["env"]
    import Pkg
    println("Activating Environment :: ", env_name)
    Pkg.activate(env_name)
end

# import packages
import NetworkInference
import LightGraphs

# load input file
if h5_flag
    println("Loading HDF5 dataset :: ", dataset_file,
            " w. ", num_genes, " genes." )
    @time genes = NetworkInference.get_h5_nodes(dataset_file,
                                                data_path=data_path,
                                                var_path=var_path,
                                                transpose_input=input_args["transpose"],
                                                number_of_nodes=num_genes);
elseif endswith(dataset_file, ".h5ad")
    println("Loading Anndata dataset :: ", dataset_file,
            " w. ", num_genes, " genes." )
    @time genes = NetworkInference.get_h5ad_nodes(dataset_file,
                                                  number_of_nodes=num_genes);
else
    println("Loading csv dataset :: ", dataset_file)
    @time genes = NetworkInference.get_nodes(dataset_file);
end
 
# generate network
println("Building Network h5ad dataset :: ", dataset_file)
algorithm = NetworkInference.PIDCNetworkInference()
@time network = NetworkInference.InferredNetwork(algorithm, genes);
 

# write output file 
println("Writing Network to file :: ", output_file)
NetworkInference.write_network_file(output_file, network);
