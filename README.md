# CeleriteQFD R package

Fitting CeleriteQFD model for detrending and flare detection in R without knowning Stan. The Stan implementation can be found [here](https://github.com/YunyiShen/CeleriteQFD). 

The model implemented here is the `CeleriteQFDexN` this model relies more on the morphology, similar to CeleriteHMM but we have three states, namely Quiet-Firing-Decay, and for **F** and **D** state, we model the data generating process as a AR(1) model with an exponential modified normal (exN) random walk for Firing and an exponential decay for Decay, this works well as pure firing model.


