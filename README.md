---
title: cmuts
emoji: 🧬
colorFrom: blue
colorTo: purple
sdk: docker
app_port: 7860
pinned: false
license: mit
---

# cmuts

Compute normalized reactivity profiles from MaP-seq experiments. Upload a reference FASTA and your sequencing reads, and get back an interactive report plus the profiles as HDF5 and CSV.

Each run also saves a `settings.json` holding every option it used. Load that file back into the form to repeat a run, or hand it to someone else with your reads.

See the [cmuts documentation](https://daslab.stanford.edu/cmuts) for details on the methods and outputs.
