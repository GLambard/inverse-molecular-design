# inverse molecular design for R: iqspr v2.3
It is devoted to the autonomous generation of novel organic compounds with target physicochemical properties initially constrained by the user. This package has the ambition to become an unavoidable tool in the innovation of novel materials and/or drugs with specific target properties. 

# Introduction

The structure of chemical species can be uniquely encoded in a single string of standard text characters called SMILES (Simplified Molecular Input Line Entry Specification). A very nice presentation of the SMILES notation can be found [here](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html). If one knows the SMILES of a chemical compound, its 
2D structure can be univoquely re-constructed. One of the aspect of the SMILES format is that it's particularly useful in the 
prediction of properties of compounds. The link that exists between the structure of a compound and its properties is generally called a QSPR (Quantitative Structure-Properties Relationship), and it has been widely used in cheminformatics for the design of new compounds. Generally, compounds structures are primarily investigated by chemists following a trial-and-error construction controlled by their existing knowledge of the chemistry and their intuition. The properties of the investigated 
compounds are then checked by direct experiments and/or driven by a QSPR analysis. In this kind of analysis, numerous descriptors can be build from the SMILES format. These descriptors can be represented a set of binary and/or continuous properties based on the existence of certain fragments in a molecule, or on the ability of its bonds to rotate for example. An introduction and overview concerning the molecular descriptors can be found [here](http://www.moleculardescriptors.eu/tutorials/tutorials.htm). Then, the descriptors are parsed as input features for a given regression model to predict output properties for a list of novel compounds. This kind of reconstruction of the properties of compounds from descriptors is called a forward prediction. 

__This package is entirely devoted to__ the inverse problem which is the backward prediction. __The generation of entirely novel SMILES in output, and consequently chemical compounds, from input targeted properties initially constrained by a user__. Thanks to the inverse-QSPR model([via](https://link.springer.com/article/10.1007%2Fs10822-016-0008-z)), this is now possible. This package has the ambition to become a useful tool, in the innovation of novel materials, in a field that is now widely referred as Materials Informatics. It is important to note that as SMILES is the basis format for this package, only organic molecular non-crystalline compounds can be generated. 

# Let's get started

* Install R >= 3.3.3 from [here](https://www.r-project.org/)

* Install RStudio from [here](https://www.rstudio.com/products/rstudio/download/#download)

* Install JAVA JDK >= 1.7 from [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)

(for issues concerning the intallation of ```rJava```, a dependency of the ```rcdk``` package included in ```iqspr```, on MAC OS X, please follow these links [here](https://github.com/snowflakedb/dplyr-snowflakedb/wiki/Configuring-R-rJava-RJDBC-on-Mac-OS-X) and [here](http://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite))

* Install the OpenBabel >= 2.3.1 with headers from [here](http://openbabel.org).

OpenBabel is compulsory to the check of the validity for the generated SMILES and their re-ordering. 

* __(Optional)__ Install ```mxnet``` >= 0.9.4 for the deep learning capabilities from [here](http://mxnet.io/get_started/setup.html).

* Install the ```devtools``` package in RStudio
```R
install.packages("devtools")
library(devtools)
```

* __(Recommended)__ Install ```rcdk``` == 3.3.8 and dependencies for better compatibilities in RStudio
```R
install_version("rcdk", version="3.3.8", repos="http://cran.us.r-project.org")
```

* Finally, install ```iqspr``` in RStudio
```R
install_github("GLambard/inverse-molecular-design",subdir="iqspr")
```

That's it! You can now use ```iqspr``` by following the tutorial delivered with the package (in *Packages* tab, click on the ```iqspr``` package in RStudio. Then, refer to the *user guide, package vignettes and other documentation*.)

# How does it work

The __iqspr__ package takes initial datasets of SMILES with their known physico-chemical properties (HOMO-LUMO gap, internal energy, melting point, toxicity, solubility, etc.) as input. Then, the SMILES are transformed in their corresponding descriptors to construct a vector of features per compound. Linear or non-linear regression models are then trained with these vectors in input, and given properties in output, to form the forward prediction model. This done, the natural language processing principle is used to build n-grams from a list of known SMILES to build a chemical grammar. Once the model for the chemical grammar is 
formed, a generator of SMILES is then available. This generator corresponds to the prior knownledge about viable chemical compounds, i.e. chemically possible, stable and which tend to be synthesizable. Finally, following the Bayes law, prior knowledge and forward prediction models, i.e. the likelihoods of a chemical structure to possess a certain property, can be linked to emulate the posterior, i.e. the probability that a given property can be represented by a given structure. Technically, thanks to a SMC (sequential Monte-Carlo), the prior distribution of possible structures is sampled (SMILES are sequentially modified character-by-character), the properties of the generated structures are predicted via the forward model, and these structures are then again transformed according to the distance of their properties to the target properties space. 
__To resume, thanks to a character-wise directed modification of SMILES, according to a prior knownledge of realistic chemical compounds, coupled to the forward model predictions, entirely novel SMILES with desired properties are autonomously generated.__
