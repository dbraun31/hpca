# Hierarchical PCA for Mind-Wandering Data

This repository contains R functions for performing hierarchical Principal
Component Analysis (PCA) on self-report data related to mind-wandering. The
methods are designed to group individuals based on the structure of their
thought patterns as inferred from self-report items.

## Overview

These functions were developed to analyze thought-probe data collected in a lab
environment. By applying hierarchical PCA, the tools identify groups of
individuals who exhibit similar patterns of thought, based on responses to
self-report items.

Key features include:

* Subject-level PCA for dimensionality reduction.
* Hierarchical clustering to group participants.
* Visualization tools for interpreting results (e.g., dendrograms, silhouette plots, and word clouds).

## Installation

Clone this repository and ensure the required R packages are installed:

```r
install.packages(c(
  "tidyr", "purrr", "stringr", "dplyr", "ggplot2", 
  "cluster", "psych", "tidygraph", "ggraph", 
  "scales", "ggwordcloud", "fs"
))
```

## Functions

### Core Functions

`HPCA(items, id = 'subject', Nfactors = 4, data = NULL, file = NULL)`

Performs hierarchical PCA on specified items.

* Parameters:

	* `items`: List of column names to include in PCA.
	* `id`: Column name for subject IDs (default: "subject").
	* `Nfactors`: Number of components to retain (default: 4).
	* `data`: Preloaded dataset (optional).
	* `file`: Path to a CSV file containing the dataset (optional).

* Returns:

	* A list containing:

		* `tree`: Dendrogram object.
		* `pcas`: Subject-level PCA results.
		* `distances`: Distance matrix.
		* `hc`: Hierarchical clustering object.

`get_ideal_groups(hpca, plot_result = TRUE, print_result = TRUE)`

Determines the optimal number of clusters using silhouette analysis.

* Parameters:

	* `hpca`: Output from HPCA().
	* `plot_result`: Display a silhouette plot (default: TRUE).
	* `print_result`: Print the ideal number of clusters (default: TRUE).

* Returns:
 
	* Ideal number of clusters.

`plot_dendro(hpca, ngroups = NULL)`

Plots a dendrogram with optional color-coded clusters.

* Parameters:

    * `hpca`: Output from HPCA().
    * `ngroups`: Number of clusters to visualize (optional).

`plot_wordcloud(hpca, ngroups, max_size = 30)`

Generates a word cloud showing PCA loadings grouped by clusters.

* Parameters:
    * `hpca`: Output from HPCA().
    * `ngroups`: Number of clusters.
    * `max_size`: Maximum text size in the word cloud (default: 30).

`write_pcas(hpca, out_path=NULL)`

Writes subject-level rotation matrices and trial-wise scores to two
separate csv's.

* Parameters:

    * `hpca`: Output from HPCA().
    * `out_path`: Directory path to save the output files. Defaults to the
        working directory.

### Auxiliary Functions

These helper functions support the main workflow and are not typically called directly:

* `rotations_to_df()`: Converts PCA loadings to a data frame.
* `get_subject_numbers()`: Extracts subject numbers from labels.
* `color_branches()`: Colors dendrogram branches by cluster.
* `pca_subjects()`: Performs PCA at the subject level.

Example Workflow

1. Load your dataset into R or provide the path to a CSV file.
2. Specify the items (columns) to include in the PCA.
3. Run HPCA() to generate the hierarchical PCA output.
4. Use get_ideal_groups() to determine the optimal number of clusters.
5. Visualize results using plot_dendro() and plot_wordcloud().

## Example

```r
# Run hierarchical PCA
results <- HPCA(items = c("item1", "item2", "item3"), file = "data.csv")

# Determine optimal clusters
ideal_clusters <- get_ideal_groups(results)

# Plot dendrogram
plot_dendro(results, ngroups = ideal_clusters)

# Generate word cloud
plot_wordcloud(results, ngroups = ideal_clusters)
```
## Notes

1. Data Normalization: Responses are normalized within each subject before PCA.
2. Clustering: The clustering algorithm uses the PCA results to compute distances between individuals.

## Contact

For questions or feedback, feel free to open an issue or reach out.
