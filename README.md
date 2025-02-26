[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/haoyunl2/SeisFlowBench/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14927938.svg)](https://doi.org/10.5281/zenodo.14927938)

# SeisFlowBench: A Benchmarking Dataset for Seismic and Flow Simulations  

This repository provides a **reproducible scientific project** built using **Julia** and **DrWatson.jl**.  
It is designed for **benchmarking seismic and flow simulations**, ensuring consistency across computational environments.  

## 📌 Project Overview  
- **Project Name:** `SeisFlowBench`  
- **Authors:** Haoyun Li and Abhinav Gahlot  
- **Focus Areas:**  
  - **Seismic Wave Propagation**  
  - **Fluid Flow in Porous Media**  
  - **Seismic-Flow Coupling Simulations**  
- **Technologies Used:**  
  - [Julia](https://julialang.org/) – A high-performance programming language for scientific computing.  
  - [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) – A framework for reproducible scientific projects.  

---

## 🚀 Getting Started  
To **reproduce this project locally**, follow these steps:  

### 1️⃣ Clone the Repository  
Download the codebase using Git:  
```sh
git clone https://github.com/your-repo/SeisFlowBench.git
cd SeisFlowBench
```  
🚨 **Note:** Raw data files are **not included** in version control (Git) and may need to be downloaded separately.  

### 2️⃣ Set Up the Julia Environment  
Open a **Julia REPL** and run the following commands:  
```julia
using Pkg
Pkg.add("DrWatson")   # Install globally (needed for `quickactivate`)
Pkg.activate(".")      # Activate this project's environment
Pkg.instantiate()      # Install all required dependencies
```  
✅ This will **install all necessary Julia packages** for the project to run correctly.  

---

## 🛠 Using the Project  
Many scripts in this repository begin with:  
```julia
using DrWatson
@quickactivate "SeisFlowBench"
```  
This ensures that:  
- The **correct project environment is activated**.  
- **Local paths are managed automatically** using DrWatson.  

After activation, you can run any script or experiment **without manually setting paths**.  

---

## 📂 Project Structure  
```
SeisFlowBench/
│── data/           # Processed and simulation data (not included in Git)
│── scripts/        # Main simulation and analysis scripts
│── src/            # Core project modules and functions
│── notebooks/      # Jupyter/Pluto notebooks for interactive analysis
│── Project.toml    # Julia package dependencies
│── Manifest.toml   # Exact dependency versions for reproducibility
│── README.md       # This file
```  

---

## 🔄 Converting JLD2 to HDF5  
Simulation results are stored in **JLD2 format**, but you may need them in **HDF5** for better compatibility.  

### 📌 Convert JLD2 to HDF5  
To convert all samples into **a single HDF5 file**, use the function provided in:  
📄 **`scripts/convert_jld2_hdf5.jl`**  

Run the script as follows:  
```julia
include("scripts/convert_jld2_hdf5.jl")
convert_all_samples_to_hdf5("2D_perm_wise", 1:10, "all_samples.h5")
```
✅ This will save all data in **one HDF5 file** under `/data/2D_perm_wise/all_samples.h5`.

---

## 📥 Loading Data from HDF5  
To **load and inspect simulation data** from HDF5, use the provided function in:  
📄 **`scripts/load_hdf5_data.jl`**  

Example usage:  
```julia
include("scripts/load_hdf5_data.jl")
hdf5_file = datadir("2D_perm_wise", "all_samples.h5")
data = load_hdf5_data(hdf5_file)

# Accessing a specific variable
S_array_sample_1 = data["sample_1"]["S_array"]
println("Loaded S_array from sample_1: Size = ", size(S_array_sample_1))
```

✅ This allows easy retrieval of seismic and flow simulation results.

---

## 📌 Additional Notes  
- If **any dependency is missing**, run `Pkg.instantiate()` again.  
- For issues, create a **GitHub issue** or contact the authors.  
- This project follows **reproducibility best practices** using DrWatson.  

---

✅ **SeisFlowBench is now ready for seismic and flow simulation benchmarking!** 🚀  
