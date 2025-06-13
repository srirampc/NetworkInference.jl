using ArgParse
using HDF5
using NetworkInference
using InformationMeasures


function write_network_csv(file_path::String, inferred_network::InferredNetwork,
                           field_sep=",")

    out_file = open(file_path, "w")

    for edge in inferred_network.edges
        nodes = edge.nodes
        write(out_file, string(
            nodes[1].label, field_sep, nodes[2].label, field_sep,
            edge.weight, "\n",
            nodes[2].label, field_sep, nodes[1].label, field_sep,
            edge.weight, "\n"
        ))
    end

    close(out_file)

end

function get_nodes_raw_data(adata::AbstractMatrix{Float32},
    gene_names::Array{String}, discretizer::String="bayesian_blocks",
    estimator::String="maximum_likelihood", number_of_bins::Int=10)
    # Initialize nodes
    number_of_nodes = size(gene_names, 1)
    nodes = Array{Node}(undef, number_of_nodes)

    for i in 1:number_of_nodes
        # nodes[i] = Node(lines[i:i, 1:end], discretizer, estimator, number_of_bins)
        label = string(gene_names[i])
        raw_values = adata[i:i, :]
        binned_values = zeros(Int, length(raw_values))
        number_of_bins = InformationMeasures.get_bin_ids!(raw_values, discretizer,
            number_of_bins, binned_values)
        probabilities = InformationMeasures.get_probabilities(estimator,
            InformationMeasures.get_frequencies_from_bin_ids(binned_values,
                number_of_bins))
        # Initialize nodes
        nodes[i] = Node(label, binned_values, number_of_bins, probabilities)
    end

    return nodes
end



function get_data_h5ad(data_file::String, data_path::String,
                       var_path::String)
    # Load data path
    adata = HDF5.h5read(data_file, data_path)
    gene_names = HDF5.h5read(data_file, var_path)
    return adata, gene_names
end


function pidc_net(nodes::Array{NetworkInference.Node}, output_file::String)
    inferred_network = InferredNetwork(PIDCNetworkInference(), nodes)
    write_network_csv(output_file, inferred_network)
end

# data_file = "/nv/hswarm1/schockalingam6/data2/CIBR/datasets/simulated/gnw2K.tsv"
# nodes = get_nodes(data_file)
# inferred_network = InferredNetwork(PIDCNetworkInference(), nodes)
# write_network_file("/nv/hswarm1/schockalingam6/data2/CIBR/output/pidc/simulated.txt", inferred_network)
function pidc_network(data_file, output_file)
    nodes = get_nodes(data_file)
    inferred_network = InferredNetwork(PIDCNetworkInference(), nodes)
    write_network_file(output_file, inferred_network)
end

function parse_commandline()
    argpv = ArgParseSettings()

    @add_arg_table argpv begin
        "--data_path", "-d"
        help = "LOOM path to matrix"
        default = "/matrix"
        "--var_gene_path", "-g"
        help = "Loom path to genes"
        default = "/row_attrs/Gene"
        "in_file"
        help = "Path to input file"
        required = true
        "out_file"
        help = "Path to output file"
        required = true
    end

    return parse_args(argpv)
end

function main()
    input_args = parse_commandline()
    println("Input args : ", input_args)
    # for (arg, val) in input_args
    #    println("   $arg : $val")
    # end
    in_loom::String = input_args["in_file"]
    data_path::String = input_args["data_path"]
    var_gene_path::String = input_args["var_gene_path"]
    out_file::String = input_args["out_file"]
    @showtime adata, gene_names = get_data_h5ad(in_loom, data_path, var_gene_path)
    @showtime adata = transpose(adata)
    println("ADATA ", eltype(adata), " ", size(adata))
    println("GENES ", eltype(gene_names), " ", size(gene_names))
    @showtime nodes = get_nodes_raw_data(adata, gene_names)
    println("NODES ", eltype(nodes), " ", size(nodes), " ", nodes[1])
    @showtime pidc_net(nodes, out_file)
end

main()
