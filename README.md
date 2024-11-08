# EvAM-MDR
Example of evolutionary accumulation modelling for multidrug resistance.

For this pipeline you'll need R with libraries `readxl` (for reading data) and `dplyr` (for data wrangling). You'll also need the `HyperTraPS-CT` source code, which can be found here https://github.com/StochasticBiology/hypertraps-ct . It has rather more dependencies, detailed in that repository. Download the `HyperTraPS-CT` code and run this code `EvAM-MDR.R` from the same directory.

This example code reads in drug resistance profiles and a phylogeny for a collection of Russian tuberculosis isolates; the dataset is from [1]. It uses `HyperTraPS-CT` functionality to cast these observations as a set of evolutionary transitions, then performs accumulation modelling with HyperTraPS to investigate the underlying evolutionary dynamics and mechanisms. You can read more here https://arxiv.org/pdf/2411.00219 .

[1] Casali, N., Nikolayevskyy, V., Balabanova, Y., Harris, S.R., Ignatyeva, O., Kontsevaya, I., Corander, J., Bryant, J., Parkhill, J., Nejentsev, S. and Horstmann, R.D., 2014. Evolution and transmission of drug-resistant tuberculosis in a Russian population. Nature genetics, 46(3), pp.279-286.
