## Abstract

from [1]:
Accurate molecular modelling of the physical and chemical behavior of highly cross-linked epoxy resins at the atomistic scale is important for the design of new property-optimized materials. However, a systematic approach to parametrizing and characterizing these systems in molecular dynamics is missing. We, therefore, present a unified scheme to derive atomic charges for amine-based epoxy resins, in agreement with the AMBER force field, based on defining reactive fragments – blocks – building the network. The approach is applicable to all stages of curing, from pure liquid, to gelation, to fully cured glass. We utilize this approach to study DGEBA/DDS epoxy systems, incorporating dynamic topology changes into atomistic Molecular Dynamics simulations of the curing reaction with 127,000 atoms. We study size effects in our simulations and predict the gel point utilizing rigorous percolation theory to accurately recover the experimental data. Furthermore, we observe excellent agreement between the estimated and the experimentally determined glass transition temperatures as a function of curing rate. Finally, we demonstrate the quality of our model by the prediction of the elastic modulus, based on uniaxial tensile tests. The presented scheme paves the way for a broadly consistent approach for modelling and characterizing all amine-based epoxy resins.

## Contact

Ana-Suncana Smith and Christian R. Wick
Friedrich-Alexander-University Erlangen-Nürnberg (FAU), Faculty of Science, Department of Physics, PULS Group, Interdisciplinary Center for Nanostructured Films (IZNF), Cauerstrasse 3, 91058, Germany 

## License

Creative Commons Attribution 4.0 International

## Context

Dataset to paper:

[1] M. Livraghi, S. Pahi, Piotr Nowakowski, D. M. Smith, C. R. Wick, A.-S. Smith, "Block chemistry for accurate modeling of epoxy resins".

All input scripts and force field files presented in [1]

## Folder structure

- amberLib

    Amber lib and frcmod files for all the Amber residues necessary to reconstruct the charge fragmentation scheme for the DGEBA epoxy resin and the amine hardeners DDS, DICY and IPTA.
    Examples containing leap input scripts to construct DGEBA molecules of length n.

- lammps

    - data

        The Lammps data files for the five DGEBA/DDS MD systems, straight out of assembling by Packmol (no minimization performed yet).

    - fixBondReact

        Pre- and post-reaction topology and map files to run the curing simulations. Compatible with Lammps Stable Dec. 2018.

    - inputScripts

        - minimizationEquilibration

            The Lammps scripts for initial minimization and equilibration. 

        - crossLinking

            Input scripts to run curing simulations and their corresponding annealing cycles.

        - tg

            Input scripts to run the annealing procedure required to collect density data for the measurement of glass transition temperature. The procedure is split into three parts due to the cluster's walltime. 
	
        - postCuringAnnealing

	        Input scripts to run the post-curing annealing protocol.

        - tensileTests

            Input script to run the uniaxial deformation.

- lammps_noCharge

    Example Lammps input/data required to simulate cross-linking in the absence of Coulomb interactions.
    
- curingReaction

    optimized structures for reactants, transition structures and products for uncatalysed concerted and step wise reactions and for the alcohol catalysed step wise reactions at wB97XD/def2-TZVP level of theory.

## Funding

This research was funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 377472739/GRK 2423/1-2019 FRASCAL