=================================================================================

  ACORNS: Absolute COncentration Robustness of Networks of Shinar-feinberg type

=================================================================================

GNU Octave (https://www.gnu.org/software/octave/index) was used to develop the functions used here.



=========
Functions
=========The function acr returns a list of species (and the Shinar-Feinberg (SF) pair associated with them) with absolute concentration robustness (ACR) in a chemical reaction network (CRN), if they exist. Once ACR in a species is found, other SF-pairs associated with the species are skipped. If no species is found or the network is not of SF-type, a message appears saying so.

The output variables 'model', 'R', 'F', and 'ACR_species' allow the user to view the following, respectively:

   - Complete network with all the species listed in the 'species' field of the structure 'model'
   - Matrix of reaction vectors of the network
   - Kinetic order matrix of the network
   - List of species with ACR

acr uses the following functions:     1. init_graph
          - OUTPUT: Creates an empty structure that represents an undirected graph. The structure has the following fields: 'vertices' and 'edges'.
          - INPUT: none     2. add_vertex
          - OUTPUT: Adds a vertex to an undirected graph. This is indicated in the 'vertices' field of the structure representing the graph.
          - INPUTS:
                    - g: a structure with fields 'vertices' and 'edges'
                    - v: a string representing the vertex     3. add_edge
          - OUTPUT: Adds an undirected edge between two vertices. The vertex connected to a vertex is indicated in the subfield 'vertex' and the label for the edge is in the subfield 'label'. Both subfields are under the field 'edges' corresponding to the vertex. The field and subfields belong to the structure representing the graph.
          - INPUTS:
                    - g: a structure with fields 'vertices' and 'edges'
                    - v1, v2: strings representing the vertices connected by an undirected edge (make sure 'add_vertex' has been used to add the vertices 'v1' and 'v2' to g)

     4. vertex_component
          - OUTPUT: Returns a vector whose entries are the component numbers where each vertex of an undirected graph belongs to. The function returns an empty value if there are no vertices in the graph.
          - INPUT: g: a structure with fields 'vertices' and 'edges'     5. extend_basis
          - OUTPUT: Basis for the reaction vectors of a CRN separated into 2 sets:
                       - B1 containing the SF-pair of reactions (or one of the two reactions, if the pair is linearly dependent) that needed to be extended to a basis for the reaction vectors
                       - B2 containing the basis vectors added to B1
          - INPUTS:
                    - SF_pair1: reaction number of the first reaction in an SF-pair
                    - SF_pair2: reaction number of the second reaction in an SF-pair
                    - R: matrix of reaction vectors of the CRN
                    - basis: basis for R
                    - basis_reaction_num: reaction numbers of the reaction vectors that form the basis for R     6. R_in_span_union
          - OUTPUT: Returns a value of 1 for the variable 'binary_decomp' if R is in the union of span(B1) and span(B2), i.e., an independent binary decomposition of R is formed. The output variables 'span_B1' and 'span_B2' allows the user to view the elements (reaction numbers) in span(B1) and span(B2), respectively.
          - INPUTS:
                    - B1: an array of the reaction numbers of an SF-pair
                    - B2: an array of the reaction numbers of the vectors added to extend B1 to a basis for R                    - R: a matrix of reaction vectors of a CRN     7. deficiency_N1
          - OUTPUT: Creates a network out of a list of reactions from another network, then returns the deficiency value of the newly-formed CRN. The output variables 'model_N1' and 'delta1' allow the user to view the following, respectively:
                       - New network formed
                       - Deficiency of the new network
          - INPUTS:
                    - model: a structure representing the CRN (see details below)
                    - span_B1: list of reaction numbers of 'model' that will be used to form 'model_N1'
        This function uses the function deficiency:
          - OUTPUT: Returns the deficiency value of a CRN through the variable 'delta'.
          - INPUT: model: a structure representing the CRN (see details below; the kinetics of the network is not needed)

     8. is_PL_RDK
          - OUTPUT: Returns a value of 1 if the power law kinetic system has reactant-determined kinetics (PL-RDK), 0 otherwise.
          - INPUT:  model: a structure representing the CRN (see details below; it is assumed that the CRN has power law kinetics)


Parts of the code come from the file model_analysis.m which is part of the ERNEST toolbox for CRN theory [7]. The computation of the number of linkage classes and strong linkage classes also utilizes functions from the same toolbox, specifically in the folders @multigraph and @umultigraph. These can all be downloaded from https://www.sissa.it/fa/altafini/papers/SoAl09/.



=====
Notes
=====

     1. Make sure all 8 functions and the folders @multigraph and @umultigraph are in the same folder/path being used as the current working directory.

     2. acr2 is the same as acr but returns each species (and the SF-pair associated with it) with ACR as they are found, if they exist. Once ACR in a species is found, other SF-pairs associated with the species are skipped. The time elapsed is also shown per species found and the end of the running time. This is perfect for networks with high rank and a lot of SF-pairs. An alphabetical list of the species is also returned at the end.

     3. acr3 is the same as acr but ACR in a species is checked for each SF-pair even if the species is already determined to have ACR considering a different SF-pair.

     4. acr4 is the same as acr2 but ACR in a species is checked for each SF-pair even if the species is already determined to have ACR considering a different SF-pair.



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

Note that for the functions acr, acr2, acr3, and acr4:

     1. It is assumed that the CRN has a positive equilibrium.
     2. It is also assumed that the CRN has power law kinetics.
     3. The CRN should have at least 2 species and 2 reactions (to form an SF-pair).
     4. Notes 2 and 3 imply that we assume the CRN is a power law kinetic system of SF-type.
     5. This code is based largely on [2] with modifications based on the ERRATUM explained in [5].



========
Examples
========

8 examples are included in this folder:

   - Example 1: Early STAT signaling network coupled with the receptor complex formation upon interferon induction [2]

   - Example 2: A network with 6 components [1]

   - Example 3: A simple two-species mass-action reaction system [4]

   - Example 4: Power law approximation of the pre-industrial carbon cycle model [3]

   - Example 5: A network with only one species [6]

   - Example 6: Another network with only one species [6]

   - Example 7: A network with two species [6]

   - Example 8: A network with power law kinetics [3]



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (21 October 2021)



==========
References
==========

   [1] Eloundou-Mbebi, J.M., KŸken, A., Omranian, N., Kleesen, S., Neigenfind, J., Basler, G., and Nikoloski, Z. (2016). A network property necessary for concentration robustness. Nature Communications, 7(13255), 1-7. doi:10.1038/ncomms13255.

   [2] Fontanil, L.L., Mendoza, E.R., and Fortun, N.T. (2021). A computational approach to concentration robustness in power law kinetic systems of Shinar-Feinberg type. MATCH Communications in Mathematical and in Computer Chemistry, 86, 489-516.

   [3] Fortun, N.T. and Mendoza, E.R. (2020). Absolute concentration robustness in power law kinetic systems. MATCH Communications in Mathematical and in Computer Chemistry, 85, 669-691.

   [4] Kuwahara, H., Umarov, R., Almasri, I., and Gao, X. (2017). ACRE: Absolute concentration robustness exploration in module-based combinatorial networks. Synthetic Biology, 2(1), 1-6. doi:10.1093/synbio/ysk01.

   [5] Lao, A.R., Lubenia, P.V.N., Magpantay, D.M., and Mendoza, E.R. (2021). Concentration robustness in LP kinetic systems (in preparation).

   [6] Meshkat, N., Shiu, A., and Torres, A. (2021). Absolute concentration robustness in networks with low-dimensional stoichiometric subspace. Vietnam Journal of Mathematics. https://doi.org/10.1007/s10013-021-00524-5.

   [7] Soranzo, N. and Altafini, C. (2009). ERNEST: a toolbox for chemical reaction network theory. Bioinformatics, 25(21), 2853Ð2854. doi:10.1093/bioinformatics/btp513.