# MCMC-MCMCABC-SMCABC

## A comparison of approximate versus exact techniques for Bayesian parameter inference in nonlinear ordinary differential equation models
### Abstract
  
The behaviour of many processes in science and engineering can be accurately described by dynamical system models consisting of a set of ordinary differential equations (ODEs). Often these models have several unknown parameters that are difficult to estimate from experimental data, in which case Bayesian inference can be a useful tool. In principle, exact Bayesian inference using Markov chain Monte Carlo (MCMC) techniques is possible; however, in practice, such methods may suffer from slow convergence and poor mixing. To address this problem, several approaches based on approximate Bayesian computation (ABC) have been introduced, including Markov chain Monte Carlo ABC (MCMC ABC) and sequential Monte Carlo ABC (SMC ABC). While the system of ODEs describes the underlying process that generates the data, the observed measurements invariably include errors. In this paper, we argue that several popular ABC approaches fail to adequately model these errors because the acceptance probability depends on the choice of the discrepancy function and the tolerance without any consideration of the error term. We observe that the so-called posterior distributions derived from such methods do not accurately reflect the epistemic uncertainties in parameter values. Moreover, we demonstrate that these methods provide minimal computational advantages over exact Bayesian methods when applied to two ODE epidemiological models with simulated data and one with real data concerning malaria transmission in Afghanistan.



In the paper their are two examples. Computer codes are provided here for all examples, with one example per folder. A description file is provided in each folder on how the code will work.
