# Chinook

Chinook is the name I gave to a set of scripts that I use to generate datasets for the purposes of [BIOL525D Bioinformatics for Evolutionary Biologists](https://ubc-biol525d.github.io/).

The datasets are based on a simulation of local adaptation in a meta-population of Chinook salmon. The simulation is implemented SLiM, in the simulation we model the evolution of a nucleotide sequence based on the chinook salmon reference genome. Check out [Simulations/](Simulations/) for more information on the simulation.

The idea of using simulated datasets is to provide a bit more transparency to students and workshop participants. Using simulations means that we have a ground truth against which to compare things that estimate using bioinformatic approaches. For example, for a *de novo* genome assembly, we can compare a particular build back the ground truth. The hope of this approach is to make bioinformatic concepts and approaches more concrete, providing a scaffold onto which students can build their knowledge and skills.

Many of the genomic resources required to implement these scripts are not hosted here due to the size limitations of GitHub repos. I've included links to specific datasets.

For anyone who has questions, suggestions or wants to contribute, please do get in touch.
______

# Simulated datasets

The pipelines in this repository make use of published methods for simulating high-throughput sequencing data. At present, we have simulated:

* Reduced representation sequences using a double digest RAD protocol using [RADinitio](http://catchenlab.life.illinois.edu/radinitio/) (Catchen et al 2020) [SimulateReads/RAD/](SimulateReads/RAD/)

* Short-read paired-end sequences modelling reads one may obtain from an Illumina MiSeq machine using [InSilicoSeq](https://github.com/HadrienG/InSilicoSeq) (ISS; Gourle et al 2018) [SimulateReads/Illumina_ISS/](SimulateReads/Illumina_ISS/)

* Long-read sequences modelled on those one may obtains using a PacBio NovaSeq machine, using [SimLoRD](https://bitbucket.org/genomeinformatics/simlord/src/master/) (Stocker et al 2016) [SimulateReads/HiFiReads/](SimulateReads/HiFiReads/)

* **Coming soon** RNA-seq

There are many programs for simulating sequencing reads. Such programs are typically used to benchmark methods. Here, we use those methods to generate data for teaching purposes.
