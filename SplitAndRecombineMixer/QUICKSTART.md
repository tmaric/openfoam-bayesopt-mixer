# Quickstart — SplitAndRecombineMixer

## Prerequisites

Build the custom OpenFOAM function objects from the repository root:

```bash
./Allwmake
```

Source the OpenFOAM environment (once per shell session):

```bash
source $HOME/OpenFOAM/OpenFOAM-v2506/etc/bashrc
```

## Run with Snakemake

```bash
cd SplitAndRecombineMixer
snakemake -j N
```

`N` is the number of CPU cores available.  It is used in two ways:

- **Snakemake** uses it to schedule independent rules concurrently.
- **OpenFOAM** (`simpleFoam`, `scalarTransportFoam`) runs in parallel with
  `mpirun -np N`, so both the flow and scalar-transport solves use all N cores.

Outputs land in `results/`:

```
results/
├── SplitAndRecombineHydro/   — mesh, flow solution, pressureDrop.csv
├── SplitAndRecombineMixing/  — scalar transport solution, mixing.csv
└── objectives.csv            — agglomerated geometry + objectives (one row)
```

## Clean

```bash
snakemake -j 1 clean
```

Removes the entire `results/` directory.
