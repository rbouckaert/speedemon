---
author: Remco Bouckaert
level: Intermediate
title: Speedemon tutorial
subtitle: Species delimitation with SNP data
beastversion: 2.7.5
---


# Background

Species delimitation based on multi species coalescent model can be done relatively quickly in a full Bayesian setting using a tree prior that collapses species in the same group when they are genetically indistinguishable.
In this tutorial, we will use a SNP alignment of East African geckos from the *Hemidactylus fasciatus* species complex (Leaché et al, 2014).

We will use a multi-species coalescent model that integrates out the gene trees (Stolz, et al, 2021) to infer species trees.

To identify species, we use a collapse model that groups and ungroups taxa during the MCMC run (Douglas & Bouckaert, 2022), depending on how genetically close the species are.


----

# Programs used in this Exercise

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees (Bouckaert et al, 2019). This tutorial uses the BEAST2 version 2.7.5.

### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Tracer

Tracer is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.7.2.

----

# Practical: Species delimitation of East African geckos

We will set up an analysis in BEAUti using the snapper package and set up a multi Yule collapse model from the speedemon package. First, we will install these packages:

> * Start BEAUti
> * Click to the `File => Manage packages` menu item.
> * Select `snapper` in the list of packages and the click `Install` button.
> * Select `speedemon` in the list of packages and the click `Install` button.
> * Close BEAUti -- it needs to restart to pick up the new packages.


## Set up in BEAUti

Before you begin, download the alignment from [here](https://raw.githubusercontent.com/BEAST2-Dev/SNAPP/master/examples/nexus/hemi129.nex).

> Start BEAUti and select the `File => Templates => snapper` item

BEAUti should change to show it uses the Fixe Tree Analysis template.

<figure>
	<a id="fig:BEAUti1"></a>
	<img style="width:45%;" src="figures/BEAUti-snapper.png" alt="">
	<img style="width:45%;" src="figures/BEAUti-snapper2.png" alt="">
	<figcaption>Figure 1: Select the snapper template, and BEAUti changes its appearance.</figcaption>
</figure>

> Select `File => Add alignment`, and choose the file `hemi129.nex` you just downloaded.

The partition panel now shows one entry per taxon, and a mapping to species. 

<figure>
	<a id="fig:BEAUti2"></a>
	<img style="width:45%;" src="figures/BEAUti-partition0.png" alt="">
	<figcaption>Figure 2: After loading an alignment, BEAUti attemtps to automatically map taxa to species.</figcaption>
</figure>

By default, species are assigned to taxa, but since we want to do species delimitation, we claim to have no knowledge about the species assignments and put every taxon in its own species. This can be done as follows:

> Click the `Guess` button. A dialog pops up where you can select `use everything` and keep the `after first` selection with the underscore as character to separate on.

<figure>
	<a id="fig:BEAUti3"></a>
	<img style="width:45%;" src="figures/BEAUti-guess.png" alt="">
	<figcaption>Figure 3: There are many ways to help in mapping taxa to species based on the taxon name. This includes specifying the mapping in a text file, where the file uses one line per taxon, and the first word is the taxon name followed by a tab (not space!) followed by the species name. Here, we use the defaults.</figcaption>
</figure>


> Click `OK`, and the species are now uniquely assigned to taxa:

<figure>
	<a id="fig:BEAUti4"></a>
	<img style="width:45%;" src="figures/BEAUti-partition1.png" alt="">
	<figcaption>Figure 4: Species assignment after using the species configuration option under the `Guess` button.</figcaption>
</figure>


> Switch to the `Model parameters` tab.

<figure>
	<a id="fig:BEAUti5"></a>
	<img style="width:45%;" src="figures/BEAUti-model-parameters.png" alt="">
	<figcaption>Figure 5: Various model parameters that can be changed.</figcaption>
</figure>

The following can be changed:
* Coalescent rate: starting values of the coalescent rates and whether they should be estimated (recommended) or stay fixed (only use when prior information is available).
* N: dimension of Chebyshef functions used, should be power of 2 plus 1. Higher values are more accurate but slower (optional, default: 33)
* non-polymorphic: Check box if there was no pre-filtering of sites to remove all constant sites. Leave unchecked if constant sites had been removed or systematically not selected (e.g. SNP data). The likelihoods will be adjusted according. (optional, default: true)
* number of sites which were not filtered to remove constant sites: Number of sites not pre-filtered.  (default =0). This setting ignored unless non-polymorphic set to TRUE (optional, default: 0)
* mutation Only At Root: Emulate the likelihood calculation of RoyChoudhury et al (2008) which assumes that mutations occur only in the ancestral (root) population (optional, default: false). Should almost certainly be left unchecked.
* use Log Likelihood Correction: use correction of log likelihood for the purpose of calculating Bayes factors for different species assignments. There is (almost) no computational cost involved for the MCMC chain, but the log likelihood might be reported as positive number with this correction since the likelihood is not a proper likelihood any more. (optional, default: true)
* use Beta Root Prior: instead of using a uniform prior for allele frequencies at the root, use a beta root prior (optional, default: false)

> Leave all values as is, and switch to the `Priors` tab

> Change the Tree prior to `Yule Skyline Collapse`

> Click the triangle next toe `Yule Skyline Collapse` to show its options

<figure>
	<a id="fig:BEAUti6"></a>
	<img style="width:45%;" src="figures/BEAUti-priors.png" alt="">
	<figcaption>Figure 6: Priors panel.</figcaption>
</figure>


The Yule Skyline Collapse is a mixture of skyline version of Yule tree prior that integrates out birth rate parameters under a gamma prior and spike distribution on internal node heights. It has the following options:

* epsilon: collapse height value below wich taxa are considered to be the same species. It may be useful to run the analysis with different values to see how sensitive the species assignments are under different settings. See section `User guide for selecting threshold ϵ` of Douglas & Bouckaert, 2022 on how to set epsilon.
* weight: mixture weight between Yule and spike density. Can be estimated (recommended).
* birth Rate Shape: Shape of the gamma prior distribution on birth rates. A value of 2 is default and is used in the literature.
* birth Rate Rate: Rate of the gamma prior distribution on birth rates. Can be estimated (recommended).
* group Count: the number of groups used, which determines the dimension of the groupSizes parameter. If less than zero (default) 10 groups will be used, unless group sizes are larger than 30 (then group count = number of taxa/30) or less than 6 (then group count = number of taxa/6 (optional, default: -1)
* equal Epochs: if equalEpochs is false, use epochs based on groups from tree intervals, otherwise use equal sized epochs that scale with the tree height (optional, default: false)
* linked Mean: use populationMean only for first epoch, and for other epochs use the posterior mean of the previous epoch (optional, default: false)


Note that at the bottom of the screen is a button to `+ Add priors`. This allows adding monophyletic constraints and other priors.

> Leave all values as is, and switch to the `MCMC` tab. Change chain length to 1 million.

<figure>
	<a id="fig:BEAUti7"></a>
	<img style="width:45%;" src="figures/BEAUti-mcmc.png" alt="">
	<figcaption>Figure 7: MCMC settings.</figcaption>
</figure>


> Save the file to XML, say `speedemon.xml`. A copy of the file can be found [here](TODO).

```
Question: There is no tab for clock models, so a relaxed clock cannot be set up in BEAUti. Any idea why this is?
```


## Running BEAST

Unfortunately, running this analysis will take about a day.
Snapper can benefit from using threads, and it depends a bit on the computer, the data


## Analysing results


----

# Useful Links

- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) (Drummond & Bouckaert, 2015)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users)

----

# Relevant References

Bouckaert R, Vaughan TG, Barido-Sottani J, Duchêne S, Fourment M, Gavryushkina A, Heled J, Jones G, Kühnert D, De Maio N, Matschiner M. BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. PLoS computational biology. 2019 Apr 8;15(4):e1006650.

Douglas J, Bouckaert R. Quantitatively defining species boundaries with more efficiency and more biological realism. Communications Biology. 2022 Jul 28;5(1):755.

Drummond AJ, Bouckaert RR. Bayesian evolutionary analysis with BEAST. Cambridge University Press; 2015.

Leaché, A. D., M. K. Fujita, V. N. Minin, and R. Bouckaert. 2014. Species delimitation using genome-wide SNP data. Systematic Biology 63:534–542.

Stoltz M, Baeumer B, Bouckaert R, Fox C, Hiscott G, Bryant D. Bayesian inference of species trees using diffusion models. Systematic Biology. 2021 Jan;70(1):145-61.