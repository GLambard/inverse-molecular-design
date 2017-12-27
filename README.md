# inverse molecular design for R: iqspr v2.4
It is devoted to the autonomous generation of novel organic compounds with target physicochemical properties initially constrained by the user. This package has the ambition to become an unavoidable tool in the innovation of novel materials and/or drugs with specific target properties. 

# Introduction

The structure of chemical species can be uniquely encoded in a single string of standard text characters called SMILES (Simplified Molecular Input Line Entry Specification). A very nice presentation of the SMILES notation can be found [here](http://www.daylight.com/dayhtml/doc/theory/theory.smiles.html). If one knows the SMILES of a chemical compound, its 
2D structure can be univoquely re-constructed. One of the aspect of the SMILES format is that it's particularly useful in the 
prediction of properties of compounds. The link that exists between the structure of a compound and its properties is generally called a QSPR (Quantitative Structure-Properties Relationship), and it has been widely used in cheminformatics for the design of new compounds. Generally, compounds structures are primarily investigated by chemists following a trial-and-error construction controlled by their existing knowledge of the chemistry and their intuition. The properties of the investigated 
compounds are then checked by direct experiments and/or driven by a QSPR analysis. In this kind of analysis, numerous descriptors can be build from the SMILES format. These descriptors can be represented a set of binary and/or continuous properties based on the existence of certain fragments in a molecule, or on the ability of its bonds to rotate for example. An introduction and overview concerning the molecular descriptors can be found [here](http://www.moleculardescriptors.eu/tutorials/tutorials.htm). Then, the descriptors are parsed as input features for a given regression model to predict output properties for a list of novel compounds. This kind of reconstruction of the properties of compounds from descriptors is called a forward prediction. 

__This package is entirely devoted to__ the inverse problem which is the backward prediction. __The generation of entirely novel SMILES in output, and consequently chemical compounds, from input targeted properties initially constrained by a user__. Thanks to the inverse-QSPR model([via](https://link.springer.com/article/10.1007%2Fs10822-016-0008-z)), this is now possible. This package has the ambition to become a useful tool, in the innovation of novel materials, in a field that is now widely referred as Materials Informatics. It is important to note that as SMILES is the basis format for this package, only organic molecular non-crystalline compounds can be generated. 

# Let's get started

## Docker image (Highly recommended)

* Install Docker (latest) for your OS from [here](https://docs.docker.com/engine/installation/)
* Go to ```Preferences``` of Docker, then click on the ```Advanced``` tab, and allocae half of your ```CPUs``` and ```Memory``` ressources to Docker. Then, close the ```Preferences``` window.  
* For MAC/Windows users, you'll need to take note of the IP address allocated to the virtual machine in which a Linux and all the necessary ressources to RStudio and iqspr are installed. For this, open a terminal window to type:
```bash 
docker-machine ls
```
In output, note the IP address indicated in the ```URL``` column such as: ```tcp://<IP_address>:<...>```

If the ```URL``` column is empty and the ```STATE``` is ```Stopped```, enter: 
```bash
docker-machine start default
```
And follow the suggestions given in output. This should take few minutes for the IP allocation.  

Then, do:
```bash
docker-machine ip default 
```
that returns the IP address you previously noted. Linux/MAC users can also just use ```localhost``` in place of this   
IP address. 

Then, make a directory where your work will be saved. For example (for Linux/MAC users): 
```bash
mkdir -p /Users/<user_name>/Documents/dockerspace
```
or (Windows users):
```bash
md c:\dockerspace
```
Then, change the permission on the created directory, such like (Linux/MAC users):
```bash
chmod 777 /Users/<user_name>/Documents/dockerspace
```
or (Windows users):
```bash
icacls "c:\dockerspace" /grant Users:F
```
You will now be able to share files between the Docker container and your machine. All the output files created or readable in the Docker container will be located in this directory or its sub-directories. 

Finally, run the container via: 
```bash
sudo docker run -d -p 8787:8787 -e ROOT=TRUE -v <working_directory>:/home/rstudio/dockerspace --name iqspr_shared lambard/iqspr
```
(Windows users do not need to use ```sudo```)
```<working_directory>``` is the directory created above. Please, respect the following naming ```/home/rstudio/<dockerspace_directory>``` for the shared directory in the Docker container. 

This should take few minutes for the downloading of the Docker image for ```iqspr```. Then, open a window on your web browser, and for the address type: 
```bash
localhost:8787
```
or, 
```bash
<IP_address>:8787
```
This should open a virtual session of RStudio in your browser. To log in, just use ```rstudio``` as username and password.

That's it! You can now use iqspr by following the tutorial delivered with the package (in Packages tab, click on the iqspr package in RStudio. Then, refer to the user guide, package vignettes and other documentation.)

### For Docker image updates
You'll have to stop and erase the running ```iqspr``` container, then erase the ```iqspr``` Docker image, to finally install the latest version, as follows:
* Note the container ID of the currently running ```iqspr``` image with:
```bash
sudo docker ps 
```
(if the container is already stopped, ```sudo docker ps -a``` will do the job)
* Then, stop and erase the container with:
```bash
sudo docker stop <container_ID>
```
Followed by, 
```bash
sudo docker rm <container_ID>
```
In this process, your working directory ```dockerspace``` is not affected and can be re-used without loss of files. 
Then, erase the current image from Docker with:
```bash
sudo docker images
```
to note the image's name, and
```bash
sudo docker rmi <image_name>
```
```<image_name>``` should be ```lambard/iqspr``` except if you tagged it with another name. 

Finally, re-run the ```iqspr``` Docker image with:
```bash
sudo docker run -d -p 8787:8787 -e ROOT=TRUE -v <working_directory>:/home/rstudio/dockerspace --name iqspr_shared lambard/iqspr
```
The latest ```iqspr``` Docker image is installed and ready. 

## Install from Source (iqspr v2.3 (obsolete) - highly depends on your system architecture)

* Install R >= 3.3.3 from [here](https://www.r-project.org/)

* Install RStudio from [here](https://www.rstudio.com/products/rstudio/download/#download)

* Install JAVA JDK <= 1.8 from [here](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)

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
