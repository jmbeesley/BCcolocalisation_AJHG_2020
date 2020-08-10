# BCcolocalisation_AJHG_2020
Scripts used for BCAC/GTEx colocalisation analysis 


#### Background 
The causal variants and genes underlying breast cancer risk associations are largely unknown. We used genetic colocalisation analysis to identify loci at which gene expression could potentially explain breast cancer phenotypes. Using data from the Breast Cancer Association Consortium and quantitative trait loci (QTL) data from the GTEx Project version 8, our findings reveal novel genetic associations between cancer phenotypes and effector genes. These results provide a source of loci for gene validation experiments.
Abundant gene expression associations mean that demonstration of genetic colocalisation is required to mitigate false positive or co-incident interpretation. We used *hyprcoloc* (C. N. Foley, J. R. Staley, P. G. Breen, B. B. Sun, P. D. W. Kirk, S. Burgess, J. M. M. Howson, A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. bioRxiv (2019), p. 592238.) to identify co-localisation between breast cancer risk (overall and estrogen receptor-negative) and gene expression in normal human breast tissue from 396 samples from GTEx v8. 


A tar.gz archive containing formatted BCAC and GTEx association data to be used as input into "run_hyprcoloc.R" can be downloaded from  https://osf.io/w7fps/download
