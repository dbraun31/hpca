source('hpca.r')

# ~~~ RUN HIERARCHICAL PCA ~~~ #

# Define items to conduct analysis on
items <- c('att', 'past', 'fut', 'self', 'ppl', 'arou', 'mvmt', 'ling', 'aff')
    
# Define file path or import the data directly
f <- 'mweeg_es.csv'
d <- read.csv('mweeg_es.csv')

# Run the analysis
# Note, groups are not currently in order of similarity
hpca <- HPCA(file='mweeg_es.csv', items=items)

# Determine ideal number of groups
groups <- get_ideal_groups(hpca)

# Visualize dendrogram and wordcloud (by group)
plot_dendro(hpca, 4)
# if you get the error 'invalid graphics state', try running: dev.off()
plot_wordcloud(hpca, 4)
