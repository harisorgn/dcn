# dcn
This is an implementation of a biophysical model of deep cerebellar nuclei projection neuron with 516 compartments, following the Hudgkey-Huxley formalism, in Julia v0.6. It was first presented in Steuber et al. 2011 and its original implementation in GENESIS is in https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=136175&amp;file=/cn_neuron_genesis_steuber_et_al/README#tabs-2.

A reduction in the complexity of the model was achieved by reducing the number of compartments following the work of Marasco et al. 2013 and implementing Strahler's analysis. The spiking dynamics were closely matched between the detailed and reduced models, integration time was greatly reduced, however the spike timings were not accurate in the reduced model as of yet.

References :

Steuber, V. et al., 2011. Determinants of synaptic integration and heterogeneity in rebound firing explored with data-driven models of deep cerebellar nucleus cells. Journal of computational neuroscience, 30(3), pp.633–658.

Marasco, A., Limongiello, A. & Migliore, M., 2013. Using Strahler’s analysis to reduce up to 200-fold the run time of realistic neuron models. Scientific reports, 3, p.2934.
