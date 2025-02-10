# EvAM-MDR
Example of evolutionary accumulation modelling for multidrug resistance.

![image](https://github.com/user-attachments/assets/fdcf8973-cc91-4ff4-a048-0a7271cf7518)

Articles
----

Here's an overview of EvAM applications to MDR, which this repo accompanies [1]. The HyperTraPS-CT machinery is described here [2].

Dependencies
----

For this pipeline you will need several R packages. The following code will install all of them.

```
## For the code in this repository
install.packages(c("readxl", "dplyr"))

## For HyperTraPS-CT itself
install.packages(c("Rcpp", "ggplot2", "ggpubr", "ggraph", "ggwordcloud", "igraph", "stringr", "stringdist", "phangorn", "phytools", "markdown"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")
```

You will also need git, wget, and a C++ compiler (Rcpp will compile the HyperTraPS-CT code). If Rcpp is working on your system then you very likely have a C++ compiler. If you don't have git, there is a workaround below (*). If you don't have wget, there is also a workaround below (**).

Running the code
----

First, clone this repository (e.g., `git clone https://github.com/StochasticBiology/EvAM-MDR.git`) or download the ZIP file and uncompress (go to the green "Code" on the upper right) it, in a directory of your choice.

Now, open the file `EvAM-MDR.R` and run it. One of the first lines in that file is a git clone command that will clone the HyperTraPS-CT repository, which can be found here https://github.com/StochasticBiology/hypertraps-ct . This cloning of the repo generally takes between 10 and 50 seconds. (*) If you don't have git, download the HyperTraPS-CT repository from that link and unzip it into a subdirectory called `hypertraps-ct` under the directory where this repo's code is sitting.

`EvAM-MDR.R` uses `wget` to download drug resistance profiles and a phylogeny for a collection of Russian tuberculosis isolates; the dataset is from [1]. (**) If you don't have wget, manually download the files from the URLs in the script and save them to the same directory as the repo's code.

`EvAM-MDR.R` uses HyperTraPS-CT functionality to cast these observations as a set of evolutionary transitions, then performs accumulation modelling with HyperTraPS to investigate the underlying evolutionary dynamics and mechanisms. You can read more here https://arxiv.org/pdf/2411.00219 . The inference code may take 1-2 hours and will output several inline messages about the state of the inference process and any numerical issues that arise during the fitting process. The remaining code analyses the output and produces several example plots.

`post-inference.Rdata` contains an R workspace image after the inference process has been run on these data; loading this should allow reproduction of the plots without re-running HyperTraPS-CT. Specifically, in file `EvAM-MDR.R`, after line 17, you can jump to line 55. Then load the above file: `load("post-inference.Rdata")`, and continue execution with line 56.

References
-----

[1] Renz, J., Dauda, K.A., Aga, O.N., Diaz-Uriarte, R., Löhr, I.H., Blomberg, B. and Johnston, I.G., 2024. Evolutionary accumulation modelling in AMR: machine learning to infer and predict evolutionary dynamics of multi-drug resistance. arXiv preprint arXiv:2411.00219. https://arxiv.org/pdf/2411.00219 
[2] Aga, O.N., Brun, M., Dauda, K.A., Diaz-Uriarte, R., Giannakis, K. and Johnston, I.G., 2024. HyperTraPS-CT: Inference and prediction for accumulation pathways with flexible data and model structures. PLOS Computational Biology, 20(9), p.e1012393. https://doi.org/10.1371/journal.pcbi.1012393
[3] Casali, N., Nikolayevskyy, V., Balabanova, Y., Harris, S.R., Ignatyeva, O., Kontsevaya, I., Corander, J., Bryant, J., Parkhill, J., Nejentsev, S. and Horstmann, R.D., 2014. Evolution and transmission of drug-resistant tuberculosis in a Russian population. Nature genetics, 46(3), pp.279-286. https://doi.org/10.1038/ng.2878
