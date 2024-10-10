# toygenerator

`toygenerator` is a Fortran Monte Carlo program designed to simulate e+ e- --> mu+ mu- (gamma) scattering events and save the results and total cross section.

## Author

- [@lorenzopunzi](https://github.com/LorenzoPunzi)

## Features
- Can generate events with (ISR) and without (Born) an ISR photon emitted from the initial state
- Can generate events in a weighted or unweighted manner
- Can apply specific cuts to kinematical quantities
- Can print the events and histograms of some quantities to output files

## Prerequisites

Before you can use `toygenerator`, ensure that you have the following installed:
- gfortran

## Run

The program can be run, after building with make, like so:

```bash
  make
  ./toygenerator path/to/input/card
```
where the path to the input card must be given.

## Help

For help on how to run the generator, use:

```bash
  make
  ./toygenerator -h
```

or 

```bash
  make
  ./toygenerator --help
```

