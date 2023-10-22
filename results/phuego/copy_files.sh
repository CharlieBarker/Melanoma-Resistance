#!/bin/bash

# List of factors and actions (increased and decreased)
factors=("Factor1" "Factor2" "Factor3")
actions=("increased" "decreased")

# Source and target directories
src_dir="results"
target_dir="fisher_summary"

# Iterate through factors and actions
for factor in "${factors[@]}"; do
    for action in "${actions[@]}"; do
        src_path="$src_dir/$factor/$action/rwr_fisher/Bfisher.txt"
        target_path="$target_dir/$factor/rwr_$action.txt"

        # Create the target directory if it doesn't exist
        mkdir -p "$(dirname "$target_path")"

        # Check if the source file exists before copying
        if [ -f "$src_path" ]; then
            cp "$src_path" "$target_path"
            echo "Copied $src_path to $target_path"
        else
            echo "Source file $src_path not found, skipping..."
        fi
    done
done

# Source and target directories
src_dir="results"
target_dir="KDE_0.5"

# Iterate through factors and actions
for factor in "${factors[@]}"; do
    for action in "${actions[@]}"; do
        src_path="$src_dir/$factor/$action/$target_dir/networks/KDE.graphml"
        target_path="$target_dir/$factor/KDE_$action.graphml"

        # Create the target directory if it doesn't exist
        mkdir -p "$(dirname "$target_path")"

        # Check if the source file exists before copying
        if [ -f "$src_path" ]; then
            cp "$src_path" "$target_path"
            echo "Copied $src_path to $target_path"
        else
            echo "Source file $src_path not found, skipping..."
        fi
    done
done

