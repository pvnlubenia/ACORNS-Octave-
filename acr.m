# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                             #
#    acr                                                                      #
#                                                                             #
#                                                                             #
# OUTPUT: Returns a list of species (and the Shinar-Feinberg (SF) pair        #
#            associated with them) with absolute concentration robustness     #
#            (ACR) in a chemical reaction network (CRN), if they exist. Once  #
#            ACR in a species is found, other SF-pairs associated with the    #
#            species are skipped. If no species is found or the network is    #
#            not of SF-type, a message appears saying so.                     #
#         The output variables 'model', 'R', 'F', and 'ACR_species' allow the #
#            user to view the following, respectively:                        #
#               - Complete network with all the species listed in the         #
#                    'species' field of the structure 'model'                 #
#               - Matrix of reaction vectors of the network                   #
#               - Kinetic order matrix of the network                         #
#               - List of species with absolute concentration robustness      #
# INPUT: model: a structure, representing the CRN, with the following fields: #
#           - id: name of the model                                           #
#           - species: a list of all species in the network; this is left     #
#                blank since incorporated into the function is a step which   #
#                compiles all species used in the model                       #
#           - reaction: a list of all reactions in the network, each with the #
#                following subfields:                                         #
#                   - id: a string representing the reaction                  #
#                   - reactant: has the following further subfields:          #
#                        - species: a list of strings representing the        #
#                             species in the reactant complex                 #
#                        - stoichiometry: a list of numbers representing the  #
#                             stoichiometric coefficient of each species in   #
#                             the reactant complex (listed in the same order  #
#                             of the species)                                 #
#                   - product: has the following further subfields:           #
#                        - species: a list of strings representing the        #
#                             species in the product complex                  #
#                        - stoichiometry: a list of numbers representing the  #
#                             stoichiometric coefficient of each species in   #
#                             the product complex (listed in the same order   #
#                             of the species)                                 #
#                   - reversible: has the value true or false indicating if   #
#                        the reaction is reversible or not, respectively      #
#                   - kinetic: has the following further subfields:           #
#                        - reactant1: a list of numbers representing the      #
#                             kinetic order of each species in the reactant   #
#                             complex in the left to right direction (listed  #
#                             in the same order of the species)               #
#                        - reactant2: a list of numbers representing the      #
#                             kinetic order of each species in the reactant   #
#                             complex in the right to left direction (listed  #
#                             in the same order of the species) (empty if the #
#                             reaction is not reversible)                     #
#        Notes:                                                               #
#           1. It is assumed that the CRN has a positive equilibrium.         #
#           2. It is also assumed that the CRN has power law kinetics.        #
#           3. The CRN should have at least 2 species and 2 reactions         #
#                 (to form an SF-pair).                                       #
#           4. Notes 2 and 3 imply that we assume the CRN is a power law      #
#                 kinetic system of SF-type.                                  #
#           5. This code is based largely on [1] with modifications based on  #
#                 the ERRATUM explained in [2].                               #
#                                                                             #
# References:                                                                 #
#   [1] Fontanil, L.L., Mendoza, E.R., and Fortun, N.T. (2021). A             #
#          computational approach to concentration robustness in power law    #
#          kinetic systems of Shinar-Feinberg type. MATCH Communications in   #
#          Mathematical and in Computer Chemistry, 86, 489-516.               #
#   [2] Lao, A.R., Lubenia, P.V.N., Magpantay, D.M., and Mendoza, E.R.        #
#          (2021). Concentration robustness in LP kinetic systems (in         #
#          preparation).                                                      #
#   [3] Soranzo, N. and Altafini, C. (2009). ERNEST: a toolbox for chemical   #
#          chemical reaction network theory. Bioinformatics, 25(21),          #
#          2853–2854. doi:10.1093/bioinformatics/btp513.                      #
#                                                                             #
# Created: 22 July 2021                                                       #
# Last Modified: 21 October 2021                                              #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



function [model, R, F, ACR_species] = acr(model)
    
    %
    % STEP 0: Add to 'model.species' all species indicated in the reactions
    %
    
    % Get all species from reactants
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            model.species{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    
    % Get species from products
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            model.species{end+1} = model.reaction(i).product(j).species;
        end
    end
    
    % Get only unique species
    model.species = unique(model.species);
    
    
    
    %
    % STEP 1: Form stoichiometric matrix N (based on [3])
    %
    
    % Count the number of species
    m = numel(model.species);
    
    % Initialize the matrix of reactant complexes
    reactant_complexes = [ ];
    
    % Initialize the matrix of product complexes
    product_complexes = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complexes(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complexes(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complexes(:, end+1) = product_complexes(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complexes(:, end+1) = reactant_complexes(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Count the total number of reactions
    r = size(N, 2);
    
    
    
    %
    % Step 2: Form kinetic order matrix F
    %
    
    % Initialize matrix F
    F = [ ];
    
    % Go through each reaction
    for i = 1:numel(model.reaction)
        
        % Case 1: The reaction is NOT reversible
        if model.reaction(i).reversible == 0
            
            % Add a row of zeros
            F(end+1, :) = zeros(1, m);
            
            % Fill out the kinetic order of all the species in the reactant
            for j = 1:numel(model.reaction(i).kinetic.reactant1)
                F(end, find(strcmp(model.reaction(i).reactant(j).species, model.species))) = model.reaction(i).kinetic.reactant1(j);
            end
        
        % Case 2: The reaction is reversible
        else
            
            % Add a row of zeros
            F(end+1, :) = zeros(1, m);
            
            % Fill out the kinetic order of all the species in the reactant in the first direction
            for j = 1:numel(model.reaction(i).kinetic.reactant1)
                F(end, find(strcmp(model.reaction(i).reactant(j).species, model.species))) = model.reaction(i).kinetic.reactant1(j);
            end
            
            % Add a row of zeros
            F(end+1, :) = zeros(1, m);
            
            % Fill out the kinetic order of all the species in the reactant in the other direction
            for j = 1:numel(model.reaction(i).kinetic.reactant2)
                F(end, find(strcmp(model.reaction(i).product(j).species, model.species))) = model.reaction(i).kinetic.reactant2(j);
            end
        end
    end
    
    
    
    %
    % STEP 3: Get the matrix of reaction vectors of the network and its rank (this point onward is based on [1])
    %
    
    % Matrix of reaction vectors
    R = N';
    
    % Rank of the network
    s = rank(N);
    
    
    
    %
    % STEP 4: Determine the linkage class where each complex/reaction belongs to
    %
    
    % Get just the unique complexes
    % ind2(i) is the index in Y of the reactant complex in reaction i
    % ind(i+r) is the index in Y of the product complex in reaction i     
    [Y, ind, ind2] = unique([reactant_complexes product_complexes]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
    
    % Initialize an undirected graph G
    G = init_graph();

    % Add vertices to G: these are the complexes
    % Note: Ci corresponds to Y(:,i)
    for i = 1:n
        
        % Use the complex number as label for each vertex
        G = add_vertex(G, strcat('C', num2str(i)));
    end
    
    % Add edges to graph G: Ci -> Cj forms an edge
    for i = 1:r
        G = add_edge(G, strcat('C', num2str([~, loc] = ismember(reactant_complexes(:,i)', Y', 'rows'))), strcat('C', num2str([~,loc] = ismember(product_complexes(:,i)', Y', 'rows'))));
    end
    
    % Determine to which linkage class (i.e., component) each complex belongs to
    linkage_class = vertex_component(G);
    
    % Initialize list of linkage class for each reaction
    reaction_linkage_class = [ ];
    
    % Label each reaction (based on its reactant complex) to which linkage class it belongs to
    for i = 1:r
        reaction_linkage_class(end+1) = linkage_class([~, loc] = ismember(reactant_complexes(:,i)', Y', 'rows'));
    end
    
    
    
    %
    % STEP 5: Get all Shinar-Feinberg pairs
    %
    
    % Initialize list of Shinar-Feinberg pairs
    SF_pair = [ ];
    SF_pair_id = { };
    
    % Go through each kinetic order vector
    for i = 1:r
        
        % Compare to each of the other kinetic order vectors
        for j = i+1:r
            
            % Get the absolute values of all entries for each row, then add the rows
            nonzeros = abs(F(i, :)) + abs(F(j, :));
            
            % Make all nonzero entries equal to 1 (for comparison later with true or false values)
            nonzeros(nonzeros ~= 0) = 1;
            
            % Check which kinetics orders match
            kinetics = F(i, :) == F(j, :);
            
            % Ignore if both kinetic orders are 0
            for k = 1:m
                if (F(i, k) == 0 && F(j, k) == 0)
                    kinetics(k) = 0;
                end
            end
            
            % Get only those rows that differ by 1 kinetic order
            if nnz(nonzeros) - 1 == sum(kinetics)
                
                % Get the index of the kinetic order that is different
                index = find(~(nonzeros == kinetics));
                
                % Compile Shinar-Feinberg pairs in matrix form: reactions 'i' and 'j' form a Shinar-Feinberg pair in species 'index')
                SF_pair(end+1, :) = [i j index];
                
                % Get the name of the species corresponding to that kinetic order
                species = model.species{index};
                
                % Compile a list of Shinar-Feinberg pairs
                SF_pair_id{end+1} = ['(R' num2str(i) ', R' num2str(j) ') in ' species];
            end
        end
    end
    
    % If there are no Shinar-Feinberg pairs, exit the algorithm
    if size(SF_pair, 1) == 0
        ACR_species = { };
        fprintf([model.id ' is not of Shinar-Feinberg type. The algorithm cannot be used. \n\n']);
        return
    end
    
    
    
    %
    % STEP 6: Form a basis for the rowspace of R
    %
    
    % Write R in reduced row echelon form: the transpose of R is used so 'basis_reaction_num' will give the pivot rows of R
    %    - 'A' is R in reduced row echelon form
    %    - 'basis_reaction_num' gives the row numbers of R which form a basis for the rowspace of R
    [A, basis_reaction_num] = rref(R');
    
    % Form the basis
    basis = R(basis_reaction_num, :);
    
    
    
    %
    % STEP 7: Initialize list of species with absolute concentration robustness (ACR)
    %
    
    % Initialize checker of species with ACR (for the loop below)
    ACR_checker = [ ];
    
    % This also serves as a control if no independent binary decomposition is found: it will remain empty
    ACR_species = { };
    
    
    
    %
    % STEP 8: Get a Shinar-Feinberg pair
    %
    
    % Go through each Shinar-Feinberg pair
    for i = 1:size(SF_pair, 1)
        
        % If ACR is already found for the species corresponding to the pair, skip the pair
        % This ensures that the list of species at the end has unique elements
        if ismember(SF_pair(i, 3), ACR_checker)
            continue
        end
        
        % Get the reaction numbers of the Shinar-Feinberg pair
        SF_pair1 = SF_pair(i, 1);
        SF_pair2 = SF_pair(i, 2);
        
        
        
    %
    % STEP 9: Extend the pair to a basis for the rowspace of R
    %
        
        % Get sets 'B1' and 'B2' from the basis formed
        [B1, B2] = extend_basis(SF_pair1, SF_pair2, R, basis, basis_reaction_num);
        
        
        
    %
    % STEP 10: Check if R is in the union of span(B1) and span(B2)
    %
        
        % Form 'span_B1' and 'span_B2'
        % 'binary_decomp' is 1 if R is in the union of span(B1) and span(B2), 0 otherwise
        [binary_decomp, span_B1, span_B2] = R_in_span_union(B1, B2, R);
        
        % Initialize list of power sets of B2 (for use in case we don't get an independent binary decomposition)
        power_set_B2 = { };
        
        % Create a list of power sets of B2
        % We exclude the empty set and the full set
        for j = 1:numel(B2)-1
            power_set_B2{j} = nchoosek(B2, j);
        end
        
        % Case 1: R is NOT in the union of span(B1) and span(B2)
        if binary_decomp == 0
            
            % NOTE: In the interest of understanding the process, I intentionally did not make the succeeding iterated loops into a function even if it will be repeated later
            
            
            
    %
    % STEP 11: Transfer some elements of B2 to B1 until an independent binary decomposition is found
    %
            
            % Go through each set in the power set
            for j = 1:numel(power_set_B2)
                
                % Go through each element of a set in the power set
                for k = 1:size(power_set_B2{j}, 1)
                    
                    % Transfer elements from B2 to B1
                    B1_new = [B1 power_set_B2{j}(k, :)];
                    B2_new = B2;
                    B2_new(ismember(B2_new, power_set_B2{j}(k, :))) = [ ];
                    
                    % Check if R is in the union of span(B1) and span(B2)
                    [binary_decomp, span_B1, span_B2] = R_in_span_union(B1_new, B2_new, R);
                    
                    % Once successful in getting an independent binary decomposition
                    if binary_decomp == 1
                        
                        
                        
    %
    % STEP 12: Get the deficiency of the network induced by span(B1)
    %
                        
                        % Create network N1 'model_N1' from span(B1) and compute its deficiency 'delta1'
                        [model_N1, delta1] = deficiency_N1(model, span_B1);
                        
                        
                        
    %
    % STEP 13: Repeat from STEP 10 until the deficiency of the network induced by span(B1) is less than or equal to 1
    %
                        
                        % Case 1: The deficiency is 0
                        if delta1 == 0
                            
                            
                            
    %
    % STEP 14: For deficiency 0, check if the induced network is a power law with reactant-determined kinetics (PL-RDK) system with the SF-pair in a linkage class
    %
                            
                            is_RDK = is_PL_RDK(model_N1);
                            if is_RDK
                                
                                % Initialize checker
                                SF_pair_same_linkage_class = [ ];
                                
                                % Check if each reaction in the SF-pair belongs to the same linkage class
                                SF_pair_same_linkage_class = reaction_linkage_class(SF_pair1) == reaction_linkage_class(SF_pair2);
                                
                                % If the SF-pair is in the same linkage class
                                if SF_pair_same_linkage_class == 1
                                    
                                    % Add the species associated with the Shinar-Feinberg pair to the checker
                                    ACR_checker(end+1) = SF_pair(i, 3);
                                    
                                    % Add the species with its associated Shinar-Feinberg pair to the list 'ACR_species'
                                    ACR_species{end+1} = [model.species{SF_pair(i, 3)}, ' [from SF-pair (R', num2str(SF_pair(i,1)), ', R' num2str(SF_pair(i,2)), ')]'];
                                    
                                    % This is created so we can exit the iterated loops
                                    flag = 1;
                                    
                                    % This breaks out of the loop for k
                                    break
                                end
                            end
                        
                        % Case 2: The deficiency is 1
                        elseif delta1 == 1
                             
                          
                          
    %
    % STEP 15: For deficiency 1, check if the induced network is a PL-RDK system (NO NEED for the condition that the SF-pair is in a linkage class)
    %
    
                            is_RDK = is_PL_RDK(model_N1);
                            if is_RDK
                                
                                % Add the species associated with the Shinar-Feinberg pair to the checker
                                ACR_checker(end+1) = SF_pair(i, 3);
                                
                                % Add the species with its associated Shinar-Feinberg pair to the list 'ACR_species'
                                ACR_species{end+1} = [model.species{SF_pair(i, 3)}, ' [from SF-pair (R', num2str(SF_pair(i,1)), ', R' num2str(SF_pair(i,2)), ')]'];
                                
                                % This is created so we can exit the iterated loops
                                flag = 1;
                                
                                % This breaks out of the loop for k
                                break
                            end
                        end
                    end
                end
                
                % This gets the signal from the loop for k
                if flag == 1
                   
                   % Reset flag
                   flag = 0;
                    
                   % This breaks out of the loop for j so we can move to the next Shinar-Feinberg pair
                   break
                end
            end
      
        % Case 2: R is the union of span(B1) and span(B2)
        % We jump to STEP 12 and continue up to STEP 15
        else
            
            % STEP 12
            [model_N1, delta1] = deficiency_N1(model, span_B1);
            
            % STEP 13
            % Case 1: The deficiency is 0
            if delta1 == 0
                
                % STEP 14
                is_RDK = is_PL_RDK(model_N1);
                if is_RDK
                    SF_pair_same_linkage_class = [ ];
                    SF_pair_same_linkage_class = reaction_linkage_class(SF_pair1) == reaction_linkage_class(SF_pair2);
                    if SF_pair_same_linkage_class == 1
                        ACR_checker(end+1) = SF_pair(i, 3);
                        ACR_species{end+1} = [model.species{SF_pair(i, 3)}, ' [from SF-pair (R', num2str(SF_pair(i,1)), ', R' num2str(SF_pair(i,2)), ')]'];
                        % No break here since we don't want to exit checking the other pairs
                    end
                end

            % Case 2: The deficiency is 1
            elseif delta1 == 1
            
                % STEP 15
                is_RDK = is_PL_RDK(model_N1);
                if is_RDK
                    ACR_checker(end+1) = SF_pair(i, 3);
                    ACR_species{end+1} = [model.species{SF_pair(i, 3)}, ' [from SF-pair (R', num2str(SF_pair(i,1)), ', R' num2str(SF_pair(i,2)), ')]'];
                    % No break here since we don't want to exit checking the other pairs
                end

            % Case 3: The deficiency is greater than 1
            else
                
                % NOTE: This repeats STEPS 11-15 where we go through each set in the power set until we get a subnetwork with deficiency less than or equal to 1
                % STEP 11
                for j = 1:numel(power_set_B2)
                    for k = 1:size(power_set_B2{j}, 1)
                        B1_new = [B1 power_set_B2{j}(k, :)];
                        B2_new = B2;
                        B2_new(ismember(B2_new, power_set_B2{j}(k, :))) = [ ];
                        [binary_decomp, span_B1, span_B2] = R_in_span_union(B1_new, B2_new, R);
                        if binary_decomp == 1
                            
                            % STEP 12
                            [model_N1, delta1] = deficiency_N1(model, span_B1);
                            
                            % STEP 13
                            if delta1 == 0
                                
                                % STEP 14
                                is_RDK = is_PL_RDK(model_N1);
                                if is_RDK
                                    SF_pair_same_linkage_class = [ ];
                                    SF_pair_same_linkage_class = reaction_linkage_class(SF_pair1) == reaction_linkage_class(SF_pair2);
                                    if SF_pair_same_linkage_class == 1
                                        ACR_checker(end+1) = SF_pair(i, 3);
                                        ACR_species{end+1} = [model.species{SF_pair(i, 3)}, ' [from SF-pair (R', num2str(SF_pair(i,1)), ', R' num2str(SF_pair(i,2)), ')]'];
                                        flag = 1;
                                        break
                                    end
                                end
                            elseif delta1 == 1
                                
                                % STEP 15
                                is_RDK = is_PL_RDK(model_N1);
                                if is_RDK
                                    ACR_checker(end+1) = SF_pair(i, 3);
                                    ACR_species{end+1} = [model.species{SF_pair(i, 3)}, ' [from SF-pair (R', num2str(SF_pair(i,1)), ', R' num2str(SF_pair(i,2)), ')]'];
                                    flag = 1;
                                    break
                                end
                            end
                        end
                    end
                    if flag == 1
                        flag = 0;
                        break
                    end
                end
            end
        end
    end
    
    % Arrange the species in alphabetical order
    ACR_species = sort(ACR_species);
    
    
    
    %
    % STEP 16: Display the results
    %
    
    % Case 1: No ACR in any species
    if numel(ACR_species) == 0
        disp(['The algorithm was not able to identify absolute concentration robustness in any species for ' model.id '.']);
        fprintf('\n\n')
    
    % Case 2: An ACR on a species exists
    else
        fprintf([model.id ' has absolute concentration robustness in the following species: \n\n']);
        fprintf('%s \n', ACR_species{:});
        fprintf('\n');
    end