# Use the following command line from deepTools to average bigwigs

bigwigAverage \
    -b {bigwig_files} \
    -o {bigwig_average} \
    -p {n_cores} \
    --verbose 


# Use the following command line to compute matrix from averaged bigwigs, on desired regions
computeMatrix scale-regions \
    -R {bed_files} \
    -S {bigwig_files} \
    -o {mat_file} \
    -p {n_cores} \
    --verbose

# Use following command line to plot profiles and heatmaps 
plotHeatmap \
    -m {mat_file} \
    -o {out_file} \
    --colorMap Blues \
    --startLabel "-500" \
    --endLabel "+500" \
    --samplesLabel "GMP Naive" "GMP KP" "Mono Naive" "Mono KP" \
    --regionsLabel {regions_label} \
    --xAxisLabel "Distance from summit" \
    --verbose 
