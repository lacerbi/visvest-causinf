# visvest-causinf (VestBMS)

Visuo-vestibular causal inference project (with Kalpana Dokka, Dora Angelakis and Weiji Ma)

## To-do list:

  1. Model recovery for implicit and explicit inference models
  2. (*Kalpana*) Learning for implicit and explicit datasets
  3. Analyze two additional subjects (implicit inference models only)
  
## Model fits list:

  1. **Explicit inference:**
     - *Probabilistic fusion* (1031): optimized, sampled, postprocessed.
     - *Bayesian, fixed-criterion probability matching [both noises]* (1041): optimized, **sampling**.
     - *Bayesian, fixed-criterion deterministic [both noises]* (1061): **none**.

  2. **Implicit inference:**
     - *Forced-fusion, Bayesian, fixed-criterion deterministic [eccentric noise]* (41): **optimizing**.
     - *Forced-fusion, Bayesian, fixed-criterion deterministic [constant noise]* (141): **optimizing**.
     - *Forced-fusion, Bayesian, fixed-criterion probability-matching [eccentric noise]* (41): **none**.
     - *Forced-fusion, Bayesian, fixed-criterion probability-matching[constant noise]* (141): **none**.

  3. **Uni-sensory localization:**
     - *Standard BDT [both noises]* (10201): optimized, sampled, postprocessed.

  4. **Joint fits:**
     - *Bayesian and fixed-criterion deterministic? [eccentric noise]*: **none**.
