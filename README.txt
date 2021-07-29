=================================================================================

  ACORNS: Absolute COncentration Robustness of Networks of Shinar-feinberg type

=================================================================================

GNU Octave (https://www.gnu.org/software/octave/index) was used to develop the functions used here.



=========
Functions
=========The function acr returns a list of species with absolute concentration robustness in a chemical reaction network (CRN), if they exist. If no such species is found or the network is not of Shinar-Feinberg type, a message appears saying so. The output variables 'model', 'R', 'F', and 'ACR_species' allow the user to view the following, respectively:

   - Complete network with all the species listed in the 'species' field of the structure 'model'
   - Matrix of reaction vectors of the network
   - Kinetic order matrix of the network
   - List of species with absolute concentration robustness

acr uses the following functions:     1. extend_basis
          - OUTPUT: Basis for the reaction vectors of a CRN separated into 2 sets:
                       - B1 containing the Shinar-Feinberg pair of reactions (or one of the two reactions, if the pair is linearly dependent) that needed to be extended to a basis for the reaction vectors
                       - B2 containing the basis vectors added to B1
          - INPUTS:
                    - SF_pair1: reaction number of the first reaction in a Shinar-Feinberg pair
                    - SF_pair2: reaction number of the second reaction in a Shinar-Feinberg pair
                    - R: Matrix of reaction vectors of the CRN
                    - basis: Basis for R
                    - basis_reaction_num: Reaction numbers of the reaction vectors that form the basis for R     2. R_in_span_union
          - OUTPUT: Returns a value of 1 for the variable 'binary_decomp' if R is in the union of span(B1) and span(B2), i.e., an independent binary decomposition of R is formed. The output variables 'span_B1' and 'span_B2' allows the user to view the elements (reaction numbers) in span(B1) and span(B2), respectively.
          - INPUTS:
                    - R: a matrix of reaction vectors of a chemical reaction network
                    - B1: an array of the reaction numbers of a Shinar-Feinberg pair
                    - B2: an array of the reaction numbers of the vectors added to extend B1 to a basis for R     3. deficiency_N1
          - OUTPUT: Creates a network out of a list of reactions from another network, then returns the deficiency value of the newly-formed CRN. The output variables 'model_N1' and 'delta1' allow the user to view the following, respectively:
                       - New network formed
                       - Deficiency of the new network
          - INPUTS:
                    - model: a structure representing the CRN (see details below)
                    - span_B1: list of reaction numbers of 'model' that will be used to form 'model_N1'
        This function uses the function deficiency:
          - OUTPUT: Returns the deficiency value of a CRN through the variable 'delta'.
          - INPUT: model: a structure representing the CRN (see details below; the kinetics of the network is not needed)

     4. is_PL_RDK
          - OUTPUT: Returns a value of 1 if the power law kinetic system has reactant-determined kinetics (PL-RDK), 0 otherwise.
          - INPUT:  model: a structure representing the CRN (see details below; it is assumed that the CRN has power law kinetics)

     5. is_min_PL_NDK
          - OUTPUT: Returns a value of 1 if the power law kinetic system which has non-reactant-determined kinetics (PL-NDK) is minimally PL-NDK, 0 otherwise.
          - INPUT:  model: a structure representing the CRN (see details below; it is assumed that the CRN has power law kinetics)

Parts of the code come from the file model_analysis.m which is part of the ERNEST toolbox for chemical chemical reaction network theory [6]. The computation of the number of linkage classes and strong linkage classes also utilizes functions from the same toolbox, specifically in the folders @multigraph and @umultigraph. These can all be downloaded from https://www.sissa.it/fa/altafini/papers/SoAl09/.



====
Note
====

Make sure all 7 functions and the folders @multigraph and @umultigraph are in the same folder/path being used as the current working directory.



=================================
How to fill out 'model' structure
=================================

'model' is the input for the function acr. It is a structure, representing the CRN, with the following fields:

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the product complex (listed in the same order of the species)
        - reversible: has the value true or false indicating if the reaction is reversible or not, respectively
        - kinetic: has the following further subfields:
             - reactant1: a list of numbers representing the kinetic order of each species in the reactant complex in the left to right direction (listed in the same order of the species)
             - reactant2: a list of numbers representing the kinetic order of each species in the reactant complex in the right to left direction (listed in the same order of the species) (empty if the reaction is not reversible)

Note that for the function acr:

     1. It is assumed that the CRN has a positive equilibrium.
     2. It is also assumed that the CRN has power law kinetics.
     3. The CRN should have at least 2 species and 2 reactions (to form a Shinar-Feinberg pair).
     4. Notes 2 and 3 imply that we assume the CRN is a power law kinetic system of Shinar-Feinberg type.



========
Examples
========

8 examples are included in this folder:

   - Example 1: Early STAT signaling network coupled with the receptor complex formation upon interferon induction [2]

   - Example 2: A network with 6 components [1]

   - Example 3: A simple two-species mass-action reaction system [4]

   - Example 4: Power law approximation of the pre-industrial carbon cycle model [3]

   - Example 5: A network with only one species [5]

   - Example 6: Another network with only one species [5]

   - Example 7: A network with two species [5]

   - Example 8: A network with power law kinetics [3]



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia



==========
References
==========

   [1] Eloundou-Mbebi, J.M., KŸken, A., Omranian, N., Kleesen, S., Neigenfind, J., Basler, G., and Nikoloski, Z. (2016). A network property necessary for concentration robustness. Nature Communications, 7(13255), 1-7. doi:10.1038/ncomms13255.

   [2] Fontanil, L.L., Mendoza, E.R., and Fortun, N.T. (2020). A computational approach to concentration robustness in power law kinetic systems of Shinar-Feinberg type. arXiv:2011.02866 [physics.chem-ph].

   [3] Fortun, N.T. and Mendoza, E.R. (2020). Absolute concentration robustness in power law kinetic systems. MATCH Communications in Mathematical and in Computer Chemistry, 85(-), 669-691.

   [4] Kuwahara, H., Umarov, R., Almasri, I., and Gao, X. (2017). ACRE: Absolute concentration robustness exploration in module-based combinatorial networks. Synthetic Biology, 2(1), 1-6. doi:10.1093/synbio/ysk01.

   [5] Meshkat, N., Shiu, A., and Torres, A. (2021). Absolute concentration robustness in networks with low-dimensional stoichiometric subspace. arXiv:2105.00109v2. 1-31.

   [6] Soranzo, N. and Altafini, C. (2009). ERNEST: a toolbox for chemical chemical reaction network theory. Bioinformatics, 25(21), 2853Ð2854. doi:10.1093/bioinformatics/btp513.