using Pkg
Pkg.activate(".")  # Activate the project environment

using JLD2, HDF5, DrWatson

"""
    convert_all_samples_to_hdf5(sim_name, samples, hdf5_filename)

Converts multiple JLD2 files into a **single** HDF5 file.

# Arguments:
- `sim_name`: Simulation name (used for directory structure)
- `samples`: A range or array of sample numbers
- `hdf5_filename`: Name of the output HDF5 file

# Output:
- A single HDF5 file containing all samples under `/sample_X` groups.
"""
function convert_all_samples_to_hdf5(sim_name, samples, hdf5_filename)
    # Define HDF5 file path
    hdf5_file = datadir(sim_name, hdf5_filename)

    # Open the HDF5 file in write mode
    h5open(hdf5_file, "w") do h5
        for sample in samples
            # Define JLD2 file path
            jld2_file = datadir(sim_name, savename(@strdict(sample), "jld2"; digits=6))

            # Check if JLD2 file exists
            if !isfile(jld2_file)
                println("⚠️ Warning: File $jld2_file not found, skipping sample $sample...")
                continue
            end

            # Load data from JLD2
            data = JLD2.load(jld2_file)

            # ✅ Ensure the sample group exists
            sample_group = create_group(h5, "sample_$sample")

            # Process and save each variable
            for (key, value) in data
                if value isa Vector{Vector{Float64}}
                    # ✅ Handle nested vectors (variable-size issue)
                    sub_group = create_group(sample_group, key)  # Create a group for this variable
                    for (i, vec) in enumerate(value)
                        sub_group["step_$i"] = vec  # Save each vector separately
                    end
                elseif value isa Vector{Matrix{Float64}}
                    # ✅ Handle Vector of Matrices (stack to 3D array if possible)
                    if all(size(m) == size(value[1]) for m in value)
                        stacked_array = cat(value...; dims=3)
                        sample_group[key] = stacked_array
                    else
                        sub_group = create_group(sample_group, key)
                        for (i, mat) in enumerate(value)
                            sub_group["step_$i"] = mat  # Save each matrix separately
                        end
                    end
                else
                    # ✅ Save scalars, vectors, and arrays directly
                    sample_group[key] = value
                end
            end

            println("✅ Stored sample $sample in $hdf5_file")
        end
    end

    println("🎉 All samples stored in one HDF5 file: $hdf5_file")
end

# Example usage: Convert samples 1 to 10 into a single HDF5 file
convert_all_samples_to_hdf5("2D_perm_perturb", 1:10, "all_samples.h5")
