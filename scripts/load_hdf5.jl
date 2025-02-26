using Pkg
Pkg.activate(".")  # Activate the project environment

using HDF5

"""
    load_hdf5_data(hdf5_filename)

Loads all available samples and their variables from the HDF5 file.

# Arguments:
- `hdf5_filename`: Path to the HDF5 file.

# Output:
- Prints available samples and their variables.
- Returns a dictionary containing all loaded data.
"""
function load_hdf5_data(hdf5_filename)
    data_dict = Dict()

    h5open(hdf5_filename, "r") do h5
        println("📂 HDF5 File: $hdf5_filename")

        # List available samples
        samples = keys(h5)
        println("🔹 Samples found:", samples)

        for sample in samples
            println("📌 Loading data for: $sample")

            # Initialize a dictionary to store the sample's data
            sample_data = Dict()

            for key in keys(h5[sample])
                value = h5[sample][key]

                if isa(value, HDF5.Group)  # Handle nested groups
                    println("  📂 Found nested dataset: $key")
                    nested_data = Dict()

                    for subkey in keys(value)
                        nested_data[subkey] = read(value[subkey])
                    end

                    sample_data[key] = nested_data
                else
                    println("  📄 Loading: $key")
                    sample_data[key] = read(value)
                end
            end

            data_dict[sample] = sample_data
        end
    end

    println("✅ Data loading complete!")
    return data_dict
end

# Example usage: Load all data from the HDF5 file
hdf5_file = datadir("2D_perm_perturb", "all_samples.h5")
data = load_hdf5_data(hdf5_file)

# Example: Access S_array from sample_1
S_array_sample_1 = data["sample_1"]["S_array"]
println("Loaded S_array from sample_1: Size = ", size(S_array_sample_1))

# Example: Access a specific step in a nested dataset
if "some_nested_vector" in keys(data["sample_1"])
    step_1_data = data["sample_1"]["some_nested_vector"]["step_1"]
    println("Loaded step_1 from some_nested_vector:", step_1_data)
end
