# Helper functions for inferring a network from a data file

"""
    get_nodes(data_file_path::String; <keyword arguments>)

Gets an array of all Nodes from a data file. It is assumed that the first
line of the file is headers (which are discarded) and the subsequent lines
each represent one node, and are of the form:

Label    data_value1  data_value2 ...

though a different delimiter may be specified.

Arguments:
* `data_file_path`: path to the data file
* `delim=false`: the file's delimiter. Leave as false if it is whitespace
* `discretizer="bayesian_blocks"`: algorithm for discretizing the data
* `estimator="maximum_likelihood"`: algorithm for estimating probabilities
* `number_of_bins=10`: will be overwritten if using "bayesian_blocks"

The "maximum_likelihood" estimator is recommended for PUC and PIDC.
"""
function get_nodes(data_file_path::String; delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks",
    estimator = "maximum_likelihood", number_of_bins = 10)

    if delim == false
        lines = readdlm(open(data_file_path); skipstart = 1)
    else
        lines = readdlm(open(data_file_path), delim; skipstart = 1)
    end
    number_of_nodes = size(lines, 1)
    nodes = Array{Node}(undef, number_of_nodes)

    for i in 1:number_of_nodes
        nodes[i] = Node(lines[i:i, 1:end], discretizer, estimator, number_of_bins)
    end

    return nodes

end


"""
    get_h5ad_nodes(data_file_path::String; <keyword arguments>)

Gets an array of all Nodes from an anndata h5ad file. It is assumed that the 
`var` data entry in anndata has a column named `gene_ids`.

Arguments:
* `data_file_path`: path to the h5ad data file
* `var_column=gene_ids`: column name in the var data frame for gene identifiers
* `discretizer="bayesian_blocks"`: algorithm for discretizing the data
* `estimator="maximum_likelihood"`: algorithm for estimating probabilities
* `number_of_bins=10`: will be overwritten if using "bayesian_blocks"

The "maximum_likelihood" estimator is recommended for PUC and PIDC.
"""

function get_h5ad_nodes(data_file_path::String; var_column::String = "gene_ids",
                        discretizer = "bayesian_blocks",
                        estimator = "maximum_likelihood",
                        number_of_bins = 10,
                        number_of_nodes = 0)
    adata = readh5ad(data_file_path)
    nsize = size(adata.X, 2)
    if number_of_nodes == 0
        number_of_nodes = nsize
    end
    if number_of_nodes > nsize
        number_of_nodes = nsize
    end
    println("Data Size: ", size(adata.X, 1), " x ", size(adata.X, 2), ".")
    println("Loading ", number_of_nodes , " genes.")
    # nodes = Array{Node}(undef, number_of_nodes)

    # for i in 1:number_of_nodes
    #     nodes[i] = Node(adata.var[:, var_column][i], adata.X[1:end, i:i],
    #                     discretizer, estimator, number_of_bins)
    # end
    nodes = pmap(i -> Node(adata.var[:, var_column][i],
                           adata.X[1:end, i:i],
                           discretizer,
                           estimator,
                           number_of_bins),
                 1:number_of_nodes)

    return nodes
end

"""
    get_h5_nodes(data_file_path::String; <keyword arguments>)

Gets an array of all Nodes from an HDF5 file. It is assumed that the 
`var` data entry in anndata has a column named `gene_ids`.

Arguments:
* `data_file_path`: path to the h5 data file
* `data_path=/X`: HDF5 dataset containing the list of gene identifiers
* `var_path=/var/gene_ids`: HDF5 dataset containing the list of gene identifiers
* `discretizer="bayesian_blocks"`: algorithm for discretizing the data
* `estimator="maximum_likelihood"`: algorithm for estimating probabilities
* `number_of_bins=10`: will be overwritten if using "bayesian_blocks"

The "maximum_likelihood" estimator is recommended for PUC and PIDC.
"""

function get_h5_nodes(data_file_path::String;
                      data_path::String = "/X",
                      var_path::String = "/var/gene_ids",
                      transpose_input = false,
                      discretizer = "bayesian_blocks",
                      estimator = "maximum_likelihood",
                      number_of_bins = 10,
                      number_of_nodes = 0)
    adata = h5read(data_file_path, data_path)
    gene_names = h5read(data_file_path, var_path)
    ngenes = size(gene_names, 1)
    if transpose_input || ngenes == size(adata, 1)
        adata = transpose(adata)
    end
    nsize = size(adata, 2)
    if number_of_nodes == 0
        number_of_nodes = nsize
    end
    if number_of_nodes > nsize
        number_of_nodes = nsize
    end
    if number_of_nodes > size(adata, 2) 
        println("Number of nodes(", number_of_nodes,
                ") > Data Size (", size(adata, 2), ")")
        return nothing
    end
    println("Data Size: ", size(adata, 1), " x ", size(adata, 2), ".")
    println("Loading ", number_of_nodes , " genes.")
    # nodes = Array{Node}(undef, number_of_nodes)

    # for i in 1:number_of_nodes
    #     nodes[i] = Node(adata.var[:, var_column][i], adata.X[1:end, i:i],
    #                     discretizer, estimator, number_of_bins)
    # end
    nodes = pmap(i -> Node(gene_names[i],
                           adata[1:end, i:i],
                           discretizer,
                           estimator,
                           number_of_bins),
                 1:number_of_nodes)

    return nodes
end



"""
    write_network_file(file_path::String, inferred_network::InferredNetwork)

Writes a network file from an InferredNetwork type. Each line of the file
will contain an edge, and since networks are assumed undirected, each edge
will be written in both directions with the same weight:

...

LabelX   LabelY  WeightXY

LabelY   LabelX  WeightXY

...

Arguments:
* `file_path`: path to the output file
* `inferred_network`: network that was inferred
"""
function write_network_file(file_path::String, inferred_network::InferredNetwork)

    out_file = open(file_path, "w")

    for edge in inferred_network.edges
        nodes = edge.nodes
        write(out_file, string(
            nodes[1].label, "\t", nodes[2].label, "\t",
            edge.weight, "\n",
            nodes[2].label, "\t", nodes[1].label, "\t",
            edge.weight, "\n"
        ))
    end

    close(out_file)

end

"""
    read_network_file(file_path::AbstractString)

Reads a network file and creates an InferredNetwork type. Assumes that the input
is such that each line contains an edge and each edge is written in both
directions with the same weight:

...

LabelX   LabelY  WeightXY

LabelY   LabelX  WeightXY

...
"""
function read_network_file(file_path::AbstractString)
    mat = readdlm(file_path)[1:2:end, :]
    edges = []
    nodes = Set()

    for i in 1:size(mat,1)
        n1_label, n2_label, weight = mat[i, :]
        n1_label = string(n1_label)
        n2_label = string(n2_label)
        n1 = Node(n1_label, [], 0, [])
        n2 = Node(n2_label, [], 0, [])
        new_edge = Edge([n1, n2], weight)
        push!(edges, new_edge)
        push!(nodes, n1_label, n2_label)
    end

    nodes = [Node(n, [], 0, []) for n in nodes]
    return InferredNetwork(nodes, edges)
end

"""
    get_adjacency_matrix(inferred_network::InferredNetwork, threshold = 1.0; <keyword arguments>)

Gets an adjacency matrix given an InferredNetwork and a threshold.

Arguments:
* `inferred_network`: network that was inferred
* `threshold=0.1`: threshold above which to keep edges in the network
* `absolute=false`: interpret threshold as an absolute confidence score

If `absolute` is false, threshold will be interpreted as the percentage of edges to keep.
"""
function get_adjacency_matrix(inferred_network::InferredNetwork, threshold = 0.1; absolute = false)

    number_of_nodes = length(inferred_network.nodes)
    adjacency_matrix = zeros(Bool, (number_of_nodes, number_of_nodes))

    labels_to_ids = Dict{String,Int}()
    ids_to_labels = Dict{Int,String}()
    i = 1
    for node in inferred_network.nodes
        labels_to_ids[node.label] = i
        ids_to_labels[i] = node.label
        i += 1
    end

    number_of_edges = absolute ?
        findfirst(x -> x.weight < threshold, inferred_network.edges) - 1 :
        Int(round(length(inferred_network.edges) * threshold))

    for edge in inferred_network.edges[1 : number_of_edges]
        node1 = labels_to_ids[edge.nodes[1].label]
        node2 = labels_to_ids[edge.nodes[2].label]
        adjacency_matrix[node1, node2] = true
        adjacency_matrix[node2, node1] = true
    end

    return adjacency_matrix, labels_to_ids, ids_to_labels

end

"""
    infer_network(data_file_path::String, inference::AbstractNetworkInference; <keyword arguments>)

Infers a network, given a data file and a network inference algorithm. It
is assumed that the first line of the file is headers (which are
discarded) and the subsequent lines each represent one node, and are of
the form:

Label    data_value1  data_value2 ...

though a different delimiter may be specified.

Arguments:
* `data_file_path`: path to the data file
* `inference`: network inference algorithm (e.g. `PIDCNetworkInference()`)
* `delim=false`: the file's delimiter. Leave as false if it is whitespace
* `discretizer="bayesian_blocks"`: algorithm for discretizing the data
* `estimator="maximum_likelihood"`: algorithm for estimating probabilities
* `number_of_bins=10`: will be overwritten if using "bayesian_blocks"
* `base=2`: base for the information measures
* `out_file_path=""`: path to output file. If empty, will not write a file

The "maximum_likelihood" estimator is recommended for PUC and PIDC.
"""
function infer_network(data_file_path::String, inference::AbstractNetworkInference; delim::Union{Char,Bool} = false,
    discretizer = "bayesian_blocks", estimator = "maximum_likelihood", number_of_bins = 10, base = 2,
    out_file_path = "")

    println("Getting nodes...")
    nodes = get_nodes(
        data_file_path,
        delim = delim,
        discretizer = discretizer,
        estimator = estimator,
        number_of_bins = number_of_bins
    )

    println("Inferring network...")
    inferred_network = InferredNetwork(inference, nodes, estimator = estimator, base = base)

    if length(out_file_path) > 1
        println("Writing network to file...")
        write_network_file(out_file_path, inferred_network)
    end

    return inferred_network

end
