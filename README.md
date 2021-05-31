# quantum-inference-protocol
In a nutshell, to estimate the amount of information required to predict the future of a stationary stochastic process using quantum computational mechanics.

---

How much memory is required to capture the dynamics of a process?
What does it take to predict the futures of the process? 

The primary motivation behind computational mechanics is to shed light on these questions through a statistical analysis of a process. The process can take the form of a sequence of observables at every time step, like the outcome of a coin flip each time the coin is flipped.

Each time step of the process has a certain past and future. One might use brute force to remember every single past and its future permutation in an attempt to predict the future. Instead of remembering every single past and its future, computational mechanics applies the equivalence relation to "merge" pasts together as long as they the futures of different pasts have the same probabilities of occurrence [1-3].

The entire dynamics of the process can then be captured in these "merged pasts", which are termed as causal states. Each causal state undergoes probabilistic transitions to the next causal state while outputting an observable. These causal states and probabilistic outputs gives rise to an edge-emitting hidden Markov model known as the epsilon machine. Allowing the epsilon machine to run and collecting the observables then reconstructs a statistically faithful process. The amount of memory required to capture the process' dynamics is given by the Shannon entropy of the stationary distribution of the causal states. This amount of memory is known as the (classical) statistical complexity.

Application of quantum mechanics to the epsilon machine allows taking advantage of the superposition principle and non-orthogonality to have a reduction of the statistical complexity [4,5], i.e. less memory is required to capture a process' dynamics and to predict the future. This is done through quantising each causal state, then computing the von Neumann entropy of the set of quantum causal states.

It was shown in Ref. [6] that it is possible to bypass the need to construct the causal states to statistically estimate the quantum statistical complexity (which we called quantum statistical memory in [6]). The code Cq_generate_concatenated was used for the results in Ref. [6]. 

Cq_generate_concatenated was primarily written in MATLAB but has been translated to Python3 though it may not be optimised (yet). 

The following article shall be cited should this code be used: https://doi.org/10.1103/PhysRevA.101.032327.

References:
[1] J. P. Crutchfield and K. Young, Physical Review Letters 63, 105 (1989).
[2] C. R. Shalizi, PhD Thesis (2001).
[3] J. P. Crutchfield, Nature Physics 8, 17 (2011).
[4] M. Gu, K. Wiesner, E. Rieper, and V. Vedral, Nature Communications 3, 762 (2012).
[5] R. Tan, D. R. Terno, J. Thompson, V. Vedral, and M. Gu, European Physical Journal Plus 129:191 (2014), 10.1140/epjp/i2014-14191-2.
[6] M. Ho, M. Gu, and Thomas J. Elliott, Physical Review A 101, 032327 (2020).
