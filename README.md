This repository contains the MATLAB files to reproduce the data and figures from [**Sparse Identification of Slow Timescale Dynamics**](https://arxiv.org/abs/2006.00940) by Jason J. Bramburger, Daniel Dylewsky, and J. Nathan Kutz (2020). Computations use the publicly available SINDy architecture found at https://faculty.washington.edu/kutz/page26/ and should be stored in a folder entitled 'Util'. Fast periods are found using the sliding-window DMD method from [**Dynamic mode decomposition for multiscale nonlinear physics**](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.063311) by Daniel Dylewsky, Molei Tao, and J. Nathan Kutz (PRE, 2020) for which the associated codes can be found at **GitHub/dylewsky/MultiRes_Discovery**. 

Data files are too large to host on GitHub and are therefore available at: https://web.uvic.ca/~jbramburger7/Publications.html

The scripts associated with this repository are as follows:

- ToyModel_SINDy.m: Continuous-time discovery of a SINDy model to fit the toy model signal. Data is loaded from toy_model_data.mat. Corresponds to work in Section II.

- ToyModel_SlowForecast.m: Discovery of a discrete-time mapping for the coarse-grained evolution of the toy model data. Data is loaded from toy_model_data.mat. Corresponds to work in Section II.

- Logistic_SlowDiscovery.m: Slow timescale discovery for the singularly perturbed logistic ODE. Corresponds to work in Section III A.

- Jupiter_SlowDiscovery.m: Slow timescale discovery for the evolution of Jupiter in its orbital plane. Data is loaded from three_body_data.mat. Corresponds to work in Section III B.

- Saturn_SlowDiscovery.m: Slow timescale discovery for the evolution of Saturn in its orbital plane. Data is loaded from three_body_data.mat. Corresponds to work in Section III B.

- Chaos_SlowDiscovery.m: Slow timescale discovery for the evolution of a signal with chaotic slow dynamics. Corresponds to work in Section III C.
