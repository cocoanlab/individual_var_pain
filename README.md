<h1> Individual Variability in Brain Representations of Pain: data and code  </h1> 
This repository contains basic data and example code for the following publication: 

<b>"Individual Variability in Brain Representations of Pain"</b>, Lada Kohoutová, Lauren Y. Atlas, Christian Büchel, Jason T. Buhle, Stephan Geuter, Marieke Jepma, Leonie Koban, Anjali Krishnan, Dong Hee Lee, Sungwoo Lee, Mathieu Roy, Scott M. Schafer, Liane Schmidt, Tor D. Wager, Choong-Wan Woo, 2022, <i>Nature Neuroscience</i>, [https://doi.org/10.1038/s41593-022-01081-x] (https://doi.org/10.1038/s41593-022-01081-x)

<h4> Abstract </h4>

Characterizing cerebral contributions to individual variability in pain processing is crucial for personalized pain medicine but has yet to be done. Here, we address this problem by identifying brain regions with high versus low inter-individual variability in their relationship with pain. We trained idiographic pain-predictive models with 13 single-trial fMRI datasets (_n_ = 404, discovery set) and quantified voxel-level importance for individualized pain prediction. With 21 regions identified as important pain predictors, we examined the inter-individual variability of local pain-predictive weights in these regions. Higher-order transmodal regions, such as ventromedial and ventrolateral prefrontal cortices, showed larger individual variability, whereas unimodal regions, such as somatomotor cortices, showed more stable pain representations across individuals. We replicated this result in an independent dataset (_n_ = 124). Overall, our study identifies cerebral sources of individual differences in pain processing, providing potential targets for personalized assessment and treatment of pain. 

<h2> Repository content </h2>

* **Example codes** (in `example_codes` folder)
	* `example_code_NN2022_individual_variability_in_pain.m` - Matlab code with an example workflow of the multivariate representational similarity analysis presented in the paper. The code uses functions from the CanlabCore toolbox (available [here](https://github.com/canlab/CanlabCore>)).
* **Region masks** (in `region_masks` folder)
	* 21 masks of the important pain-predictive regions identified in the study (in NIfTI format)
* **Individualised predictive maps** (in `predictive_maps` folder)
	* Pain-predictive maps for 404 individuals from the study (in NIfTI format)
* **Normalised representational dissimilarity matrices** (in `normalised_RDMs` folder)
	*  Vectorised lower triangles of the normalised representational dissimilarity matrices for each of the 21 important pain-predictive regions (in .mat format). The figures in the paper can be reproduced using these data.





