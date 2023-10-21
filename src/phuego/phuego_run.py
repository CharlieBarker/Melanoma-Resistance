import os
from phuego import phuego

# User input: paths.
support_data_folder = "/Users/charliebarker/Desktop/phuego_support"
res_folder = "/Users/charliebarker/Desktop/Melanoma_Resistance/results/phuego"

# Define the test directory.
test_dir = "/Users/charliebarker/Desktop/Melanoma_Resistance/results/mofa/mofa_factors"

# When calling the functions, set your parameters here.
fisher_geneset = ["B"]
fisher_threshold = 0.05
fisher_background = "intact"
ini_pos = ["False"]
ini_neg = ["False"]
damping = 0.85
rwr_threshold = 0.05
kde_cutoff = [0.85]
use_existing_rwr = False
convert2folder = True
include_isolated_egos_in_KDE_net = False
net_format = "graphml"

# Create a results subfolder inside the main results folder.
results_subfolder = os.path.join(res_folder, "results")

# Ensure the results subfolder exists.
os.makedirs(results_subfolder, exist_ok=True)

# Iterate over files in the test directory.
for filename in os.listdir(test_dir):
    if filename.endswith(".txt"):
        # Remove the ".txt" extension from the filename.
        filename_without_extension = os.path.splitext(filename)[0]

        # Construct the full path to the test file.
        test_file = os.path.join(test_dir, filename)

        # Construct the result folder for this iteration.
        result_folder = os.path.join(results_subfolder, filename_without_extension)

        # Ensure the result folder exists.
        os.makedirs(result_folder, exist_ok=True)

        # Run phuego for this file.
        print(f"Running phuego on file: {filename}")
        phuego(
            support_data_folder=support_data_folder,
            res_folder=result_folder,
            test_path=test_file,
            fisher_geneset=fisher_geneset,
            fisher_threshold=fisher_threshold,
            fisher_background=fisher_background,
            ini_pos=ini_pos,
            ini_neg=ini_neg,
            damping=damping,
            rwr_threshold=rwr_threshold,
            kde_cutoff=kde_cutoff,
            use_existing_rwr=use_existing_rwr,
            convert2folder=convert2folder,
            include_isolated_egos_in_KDE_net=include_isolated_egos_in_KDE_net,
            net_format=net_format,
        )
