# Simulating sequence data and metadata for epidemic spread between 3 demes

This workflow simulates a simple sequence alignment (FASTA file) and a metadata file (csv file) from an SEIR epidemic spreading in a structured population composed of 3 subgroups (A, B, and C).
The metadata file describes the sequence label, sampling time and the subgroup of each sampled individual. 

The simulations outputs can be found at:
- [`results/simulated_alignment.fasta`](https://github.com/blab/ai-phylogeo/blob/main/simulations/results/simulated_alignment.fasta) for the simulated sequence alignment (73 sequences)
- [`results/simulated_metadata.csv`](https://github.com/blab/ai-phylogeo/blob/main/simulations/results/simulated_metadata.csv) for the corresponding simulated metadata (73 sequences)

We now have two larger simulated datasets that can be found at:
- [`results/simulated_alignment_larger.fasta`](https://github.com/blab/ai-phylogeo/blob/main/simulations/results/simulated_alignment_larger.fasta) for the simulated sequence alignment (2434 sequences)
- [`results/simulated_metadata_larger.csv`](https://github.com/blab/ai-phylogeo/blob/main/simulations/results/simulated_metadata_larger.csv) for the corresponding simulated metadata (2434 sequences)

- [`results/simulated_alignment_much_larger.fasta`](https://github.com/blab/ai-phylogeo/blob/main/simulations/results/simulated_alignment_much_larger.fasta) for the simulated sequence alignment (23095 sequences)
- [`results/simulated_metadata_much_larger.csv`](https://github.com/blab/ai-phylogeo/blob/main/simulations/results/simulated_metadata_much_larger.csv) for the corresponding simulated metadata (23095 sequences)

I have compressed those files which can be decompressed locally using the [zstd package](https://www.mankier.com/1/zstd). 

**Note:** This workflow relies on a local installation of BEAST and R. Next steps could include including them in a Docker container.