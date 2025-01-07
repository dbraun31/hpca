library(tidyr)
library(purrr)
library(stringr)
library(dplyr)
library(ggplot2)
library(cluster)
library(psych)
library(tidygraph)
library(ggraph)
library(scales)
library(ggwordcloud)

# --- CORE FUNCTIONS --- #
# use these functions to initiate and see the results of the hierarchical PCA analysis

HPCA <- function(items, id='subject', Nfactors=4, data=NULL, file=NULL) {
    # Perform hierarchical PCA on specified items using subject-level PCA and clustering.
    # 
    # Parameters
    # ----------
    # items : list of str
    #     Column names of items to include in PCA.
    # id : str, optional
    #     Column name for subject IDs. Default is 'subject'.
    # Nfactors : int, optional
    #     Number of principal components to retain. Default is 4.
    # data : data.frame, optional
    #     Pre-loaded dataset. If not provided, `file` must be specified.
    # file : str, optional
    #     Path to CSV file containing the dataset.
    # 
    # Returns
    # -------
    #     dict
    #         Dictionary containing the hierarchical clustering tree, subject-level PCAs, distance matrix, and hierarchical clustering object.
    # 
    # Raises
    # ------
    # Error
    #   If both `data` and `file` are None.
    
    # Read in data from file if no data is provided
    if (is.null(data) & is.null(file)) stop('Need to specify data or input path (file).')
    if (is.null(data)) data <- read.csv(file)
    d <- data
    rm(data)
    
    # Keep only relevant data
    exact_match <- colnames(d) %in% items
    response_match <- colnames(d) %in% paste0(items, '_response')
    id_match <- colnames(d) == id
    col_match <- exact_match | response_match | id_match
    d <- d[, col_match]
    d <- drop_suffix(d)
    
    # Normalize within subject
    d <- normalize_subject_item(d)
    
    # Split by subject
    d_split <- split(d, paste0('subject', d$subject))
    
    # Get subject-level PCAs
    pcas <- lapply(d_split, pca_subjects, Nfactors=Nfactors)
    pcas <- Filter(negate(is.null), pcas)
    subject_mapping <- get_subject_numbers(names(pcas))
    
    # Compute distances
    # rows and columns in this matrix correspond to index in pcas list
    dist_matrix <- get_distances(pcas)
    
    # Cluster
    hc <- hclust(as.dist(dist_matrix))
    dendro <- as.dendrogram(hc)
    # Relabel leaves
    dendro <- dendrapply(dendro, relabel_leaves, labels = subject_mapping)
    
    out <- list(tree = dendro, pcas = pcas, distances = dist_matrix, hc = hc)
    
    return(out)
    
}


get_ideal_groups <- function(hpca, plot_result=TRUE, print_result=TRUE) {
    # Determine the optimal number of clusters using silhouette analysis.
    # 
    # Parameters
    # ----------
    # hcpa : list
    #     Output from the HPCA function.
    # plot_result : bool, optional
    #     Whether to display a silhouette plot. Default is True.
    # print_result : bool, optional
    #     Whether to print the ideal number of clusters. Default is True.
    # 
    # Returns
    # -------
    # int
    #     Ideal number of clusters based on maximum average silhouette width.
    
    hc <- hpca[['hc']]
    dist_matrix <- hpca[['distances']]
    sil_width <- c()
    for (i in 2:(nrow(dist_matrix) - 1)) {
      cluster <- cutree(hc, k=i)
      sil <- as.data.frame(cluster::silhouette(cluster, dist(dist_matrix)))
      sil_width <- c(sil_width, mean(sil$sil_width))
    }
    
    sil_result <- data.frame(cuts = 2:(nrow(dist_matrix) - 1), sils = sil_width)
    sil_result$label <- ifelse(sil_result$sils > quantile(sil_result$sils, .80),
                               paste0('(', sil_result$cuts, ', ', round(sil_result$sils, 3), ')'),
                               '')
    
    max_cuts <- sil_result[sil_result$sils == max(sil_result$sils),]$cuts
    
    p <- sil_result %>% 
        ggplot(aes(x = cuts, y = sils)) + 
        geom_vline(xintercept = max_cuts, linetype='dashed', color = 'steelblue', size = 2) + 
        geom_line() +
        annotate('text', x = 10, y = .2, label = paste0('Ideal number of groups: ', max_cuts)) + 
        labs(
            x = 'Number of clusters',
            y = 'Average silhouette width'
        ) + 
        theme_bw() + 
        theme(
            axis.ticks = element_blank(),
            panel.grid = element_blank()
        )
    
    if (print_result) print(paste0('Ideal number of groups is ', max_cuts))
    if (plot_result) print(p)
    return(max_cuts)
    
}

plot_dendro <- function(hpca, ngroups=NULL) {
    # Plot the dendrogram from hierarchical PCA.
    # 
    # Parameters
    # ----------
    # hpca : list
    #     Output from the HPCA function.
    # ngroups : int, optional
    #     Number of clusters to visualize with colored branches. Default is None.
    # 
    # Returns
    # -------
    # None
    #     Displays the dendrogram plot.
    
    tree <- hpca[['tree']]
    hc <- hpca[['hc']]
    pcas <- hpca[['pcas']]
    
    
    if (!is.null(ngroups)) { 
        cluster_groups <- cutree(hc, ngroups)
        cluster_colors <- scales::brewer_pal(palette='Dark2')(max(cluster_groups))
        cluster_map <- setNames(cluster_colors[cluster_groups], names(pcas))
        tree <- dendrapply(tree, color_branches, cluster_map)
    }
    
    
    plot(tree, xlab = 'Subject ID')
    if (!is.null(ngroups)) {
        unique_clusters <- unique(cluster_groups)
        legend('topright',
               legend = paste('Group', unique_clusters),
               col = cluster_colors[unique_clusters],
               lwd = 2,
               title = 'Clusters')
    }
    
}

plot_wordcloud <- function(hpca, ngroups, max_size = 30) {
    # Generate a word cloud showing PCA loadings grouped by clusters.
    # 
    # Parameters
    # ----------
    # hpca : list
    #     Output from the HPCA function.
    # ngroups : int
    #     Number of clusters to display.
    # max_size : int, optional
    #     Maximum text size in the word cloud. Default is 30.
    # 
    # Returns
    # -------
    # matplotlib.figure.Figure
    #     Word cloud plot showing loadings by PCA components and groups.
        
    pcas <- hpca[['pcas']]
    hc <- hpca[['hc']]
    
    gold <- '#DDA812'
    navy <- '#10294B'
    
    groups <- cutree(hc, ngroups)
    to_merge <- data.frame(subject = unname(get_subject_numbers(names(pcas))), group = groups)
    
    d <- rotations_to_df(pcas)
    
    out <- d %>% 
        inner_join(to_merge) %>% 
        gather(PC, loading, PC1:(ncol(d))) %>% 
        group_by(group, PC, item) %>% 
        summarize(loading = mean(loading), N = n()) %>% 
        mutate(group = paste0('Group ', group)) %>% 
        ggplot(aes(label=item, size = abs(loading), color = loading)) +
        geom_text_wordcloud(area_corr = TRUE) +
        scale_size_area(max_size = max_size) +
        scale_color_gradient(low=navy, high=gold) + 
        labs(
            caption = 'Gold and navy indicate positive and negative loadings, respectively.'
        ) + 
        facet_grid(group~PC) + 
        theme_bw() + 
        theme(strip.background = element_rect(fill=NA))
    
    return(out)
}


# ~~ AUX FUNCTIONS ~~ #
# (these functions shouldn't need to be called directly)

rotations_to_df <- function(pcas) {
    # Convert pcas rotations to a data frame
    
    parse_pca <- function(subject_str) {
        ds <- pcas[[subject_str]][['rotations']]
        items <- rownames(ds)
        ds <- data.frame(subject=unname(get_subject_numbers(subject_str)), item = items, ds)
        return(ds)
    }
    
    d <- do.call(rbind, lapply(names(pcas), parse_pca))
    row.names(d) <- 1:(nrow(d))
    
    return(d)
}


get_subject_numbers <- function(labels) {
    # Get subject numbers (as int) from string labels (eg, 'subject1')

    match_fn <- function(name) as.integer(str_match(name, pattern = 'subject(\\d+)')[,2])
    out <- sapply(labels, match_fn)
    return(out)
    
}

color_branches <- function(node, cluster_map) {
    # Color branches of dendrogram by group
    
    if (is.leaf(node)) {
        node_label <- attr(node, 'label')
        node_color <- cluster_map[paste0('subject', node_label)]
        attr(node, 'edgePar') <- list(col = node_color, lwd = 2)
    }
    return(node)
}


pca_subjects <- function(subd, Nfactors) {
    # See here for how to actually derive a varimax rotation matrix
    # https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
    # Using the psych package below
    
    subject <- unique(subd$subject)
    pca_data <- subd[, ! colnames(subd) %in% c('subject', 'subject_idx')]
    if (0 %in% apply(pca_data, MARGIN=2, FUN=function(x) var(x))) {
        return(NULL)
    }
    
    pca_result <- psych::principal(pca_data, nfactors = Nfactors, rotate='varimax')
    loadings <- matrix(pca_result$loadings, ncol=Nfactors)
    rownames(loadings) <- rownames(pca_result$loadings)
    colnames(loadings) <- paste0('PC', 1:(ncol(loadings)))
    eigens <- pca_result$Vaccounted['Proportion Explained',]
    out <- list(
        'rotations' = loadings,
        'eigens' = eigens
    )
    return(out)
}


get_distances <- function(pcas) {
    
    # Construct an N x N matrix of zeros
    distance_matrix <- matrix(0, nrow=length(pcas), ncol=length(pcas))
    
    for (i in 1:nrow(distance_matrix)) {
        for (j in 1:ncol(distance_matrix)) {
            # Mean Euclidean distance between each column of two matrices weighted by eigenvalues
            # This parses through pcas by *index*, not by subject number
            distance_matrix[i,j] <- weighted_distance(pcas[[i]][['rotations']] - pcas[[j]][['rotations']],
                                                      e1=pcas[[i]][['eigens']], e2=pcas[[j]][['eigens']])
        }
    }
    
    return(distance_matrix)
}


relabel_leaves <- function(dendro, labels) {
    if (is.leaf(dendro)) {
        idx <- as.numeric(attr(dendro, 'label'))
        attr(dendro, 'label') <- labels[idx]
    }
    return(dendro)
}


weighted_distance <- function(diff_mat, e1, e2) {
    # Takes in a matrix of element-wise differences between two matrices
    # Also takes in two vectors of eigen values for the PCs of the two matrices
    # Returns a scalar representing the mean of the column-wise Euclidean distances
    # Weighted by the eigenvalues
    
    euclids <- sqrt(colSums(diff_mat^2))
    rel_var <- (e1 + e2) / sum(e1, e2)
    return(sum(euclids * rel_var))
}

normalize_subject_item <- function(d) {
    # Normalize within subject and item
    gather_cols <- colnames(d)
    gather_cols <- gather_cols[!gather_cols %in% c('subject', 'item', 'conf')]
    d <- d %>% 
        mutate(id = 1:(nrow(d))) %>% 
        gather(item, response, gather_cols[1]:gather_cols[length(gather_cols)]) %>% 
        group_by(subject, item) %>% 
        mutate(response_m = mean(response), response_sd = sd(response)) %>% 
        ungroup() %>% 
        mutate(response_sc = custom_scale(response, response_m, response_sd)) %>% 
        select(-response, -response_m, -response_sd) %>% 
        rename(response = response_sc) %>% 
        spread(item, response) %>% 
        select(-id) 
    return(d)
}

drop_suffix <- function(d) {
    # Drop everything after underscore from variable name
    out <- sapply(colnames(d), FUN=function(x) unlist(strsplit(x, '_'))[1])
    colnames(d) <- out
    return(d)
}


custom_scale <- function(response, response_m, response_sd) {
    # Check to ensure there's a non zero sd
    scaled <- (response - response_m) / response_sd
    out <- ifelse(response_sd != 0 & !is.na(response_sd), scaled, response_m)
    return(out)
}




