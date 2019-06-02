# dcn
This is an implementation of a biophysical model of deep cerebellar nuclei projection neuron with 516 compartments, following the Hudgkey-Huxley formalism, in Julia v0.6. It was first presented in Steuber et al. 2011 and its original implementation in GENESIS is in https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=136175&amp;file=/cn_neuron_genesis_steuber_et_al/README#tabs-2. 
An initial attempt was made to reduce the complexity of the model by reducing the number of compartments following the work of Marasco et al. 2013. The spiking dynamics were closely matched between the detailed and reduced models, however the spike timings are not accurate in the reduced model as of yet. 