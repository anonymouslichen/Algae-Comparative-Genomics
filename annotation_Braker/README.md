Files in this folder are scripts used to annotate 33 algal genomes using BRAKER2. Run in the following order:

repeatmasker.sh - this masks repeats, an essential step for downstream annotation with BRAKER

run_Braker.sh - this runs the annotation 

add_prefix.sh - this adds a prefix to the Braker output file names and headers within those files. These can them be used in Orthofinder to find single copy orthologous groups
