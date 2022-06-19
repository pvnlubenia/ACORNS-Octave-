% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                             %
%    acr                                                                      %
%                                                                             %
%                                                                             %
% OUTPUT: Returns a list of species (together with the Shinar-Feinberg (SF)   %
%            pair associated with each and the deficiency of the building     %
%            block subnetwork containing the SF-pair) with absolute           %
%            concentration robustness (ACR) in a chemical reaction network    %
%            (CRN), if they exist. ACR in a species is checked for each       %
%            SF-pair even if the species is already determined to have ACR    %
%            considering a different SF-pair. If no species is found or the   %
%            network is not of SF-type, a message appears saying so.          %
%         The output variables 'model', 'R', 'F', and 'ACR_species' allow the %
%            user to view the following, respectively:                        %
%               - Complete network with all the species listed in the         %
%                    'species' field of the structure 'model'                 %
%               - Matrix of reaction vectors of the network                   %
%               - Kinetic order matrix of the network                         %
%               - List of species with absolute concentration robustness      %
% INPUT: model: a structure, representing the CRN, with the following fields: %
%           - id: name of the model                                           %
%           - species: a list of all species in the network; this is left     %
%                blank since incorporated into the function is a step which   %
%                compiles all species used in the model                       %
%           - reaction: a list of all reactions in the network, each with the %
%                following subfields:                                         %
%                   - id: a string representing the reaction                  %
%                   - reactant: has the following further subfields:          %
%                        - species: a list of strings representing the        %
%                             species in the reactant complex                 %
%                        - stoichiometry: a list of numbers representing the  %
%                             stoichiometric coefficient of each species in   %
%                             the reactant complex (listed in the same order  %
%                             of the species)                                 %
%                   - product: has the following further subfields:           %
%                        - species: a list of strings representing the        %
%                             species in the product complex                  %
%                        - stoichiometry: a list of numbers representing the  %
%                             stoichiometric coefficient of each species in   %
%                             the product complex (listed in the same order   %
%                             of the species)                                 %
%                   - reversible: has the value true or false indicating if   %
%                        the reaction is reversible or not, respectively      %
%                   - kinetic: has the following further subfields:           %
%                        - reactant1: a list of numbers representing the      %
%                             kinetic order of each species in the reactant   %
%                             complex in the left to right direction (listed  %
%                             in the same order of the species)               %
%                        - reactant2: a list of numbers representing the      %
%                             kinetic order of each species in the reactant   %
%                             complex in the right to left direction (listed  %
%                             in the same order of the species) (empty if the %
%                             reaction is not reversible)                     %
%        Notes:                                                               %
%           1. It is assumed that the CRN has a positive equilibrium.         %
%           2. It is also assumed that the CRN has power law kinetics.        %
%           3. The CRN should have at least 2 species and 2 reactions         %
%                 (to form an SF-pair).                                       %
%           4. Notes 2 and 3 imply that we assume the CRN is a power law      %
%                 kinetic system of SF-type.                                  %
%           5. This code is based largely on [1] with modifications based on  %
%                 [2].                                                        %
%                                                                             %
% References                                                                  %
%    [1] Fontanil L, Mendoza E, Fortun N (2021) A computational approach to   %
%           concentration robustness in power law kinetic systems of          %
%           Shinar-Feinberg type. MATCH Commun Math Comput Chem               %
%           86(3):489-516.                                                    %
%    [2] Lao A, Lubenia P, Magpantay D, Mendoza E (2022) Concentration        %
%           robustness in LP kinetic systems. MATCH Commun Math Comput Chem   %
%           88(1):29-66. https://doi.org/10.46793/match.88-1.029L             %
%    [3] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction %
%           network theory. Bioinform 25(21):2853–2854.                       %
%           https://doi.org/10.1093/bioinformatics/btp513                     %
%                                                                             %
% Created: 22 July 2021                                                       %
% Last Modified: 19 June 2022                                                 %
%                                                                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [model, R, F, ACR_species] = acr(model)
    
    %
    % STEP 1: Add to 'model.species' all species indicated in the reactions
    %
    
    [model, m] = model_species(model);
    
    
    
    %
    % STEP 2: Form stoichiometric matrix N
    %
    
    [N, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    
    
    %
    % Step 3: Form kinetic order matrix F
    %    Reminder: Algorithm assumes the system is mass action
    %
    
    F = kin_ord_matrix(model, m);
    
    
    
    %
    % STEP 4: Get the matrix of reaction vectors of the network and its rank
    %
    
    R = N';
    
    
    
    %
    % STEP 5: Create a list of reactant complexes
    %
    
    % Get just the unique complexes
    [Y, ~, ~] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
    
    % Initialize list of reactant complexes
    reactant_complex_list = { };
    
    % Go through each reactant complex
    for i = 1:size(reactant_complex, 2)
        
        % For the zero complex
        if numel(find(reactant_complex(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(reactant_complex(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if reactant_complex(:, i)(find(reactant_complex(:, i))(j)) == 1
                        complex = [model.species{find(reactant_complex(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(reactant_complex(:, i)(find(reactant_complex(:, i))(j))), model.species{find(reactant_complex(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if reactant_complex(:, i)(find(reactant_complex(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(reactant_complex(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(reactant_complex(:, i)(find(reactant_complex(:, i))(j))), model.species{find(reactant_complex(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add the complex in the list
        reactant_complex_list{end+1} = complex;
    end
    
    
    
    %
    % STEP 6: Get all Shinar-Feinberg pairs
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
                
                % Compile Shinar-Feinberg pairs in matrix form: reactions i and j form a Shinar-Feinberg pair in species 'index')
                SF_pair(end+1, :) = [i, j, index];
                
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
    % STEP 7: Form a basis for the rowspace of R
    %
    
    % Write R in reduced row echelon form: the transpose of R is used so basis_reac_num will give the pivot rows of R
    %    - basis_reac_num gives the row numbers of R which form a basis for the rowspace of R
    [~, basis_reac_num] = rref(R');
    
    % Form the basis
    basis = R(basis_reac_num, :);
    
    % Initialize list of species with absolute concentration robustness (ACR)
    % This also serves as a control if no independent binary decomposition is found: it will remain empty
    ACR_species = { };
    
    
    
    %
    % STEP 8: Get a Shinar-Feinberg pair
    %
    
    % Initialize flag for the loops
    flag = 0;
    
    % Go through each Shinar-Feinberg pair
    for i = 1:size(SF_pair, 1)
        
        % Get the reaction numbers of the Shinar-Feinberg pair
        SF_pair1 = SF_pair(i, 1);
        SF_pair2 = SF_pair(i, 2);
        
        
        
    %
    % STEP 9: Extend the pair to a basis for the rowspace of R
    %
        
        % Note: A vector v is a linear combination of vectors in a matrix A if the rank of the matrix [A, v] (v is appended to A) is the same as the rank of A
        
        % Initialize the extended basis
        basis_SF = [ ];
        
        % Case 1: The Shinar-Feinberg pair is NOT linearly independent, i.e., each is a linear combo of the other
        if rank([R(SF_pair1, :)', R(SF_pair2, :)']) == rank(R(SF_pair1, :)')
        
            % Use the first reaction vector of the Shinar-Feinberg pair to start the formation of the basis
            basis_SF(end+1, :) = R(SF_pair1, :);
            
            % This forms B1
            B1 = SF_pair1;
            
        % Case 2: The Shinar-Feinberg pair is linearly independent
        else
        
            % Use the reaction vectors to start the formation of the basis
            basis_SF(end+1, :) = R(SF_pair1, :);
            basis_SF(end+1, :) = R(SF_pair2, :);
            
            % These form B1
            B1 = [SF_pair1, SF_pair2];
        end
        
        % Initialize list of reaction vectors added to extend basis_SF to a basis for R
        added_reaction = [ ];
        
        % Keep on going through each vector in 'basis' until there are s [= rank(N)] elements in basis_SF
        while size(basis_SF, 1) < size(basis, 1)
        
            % Go through each vector in 'basis'
            for j = 1:size(basis, 1)
            
                % If the reaction vector is NOT a linear combination of basis_SF
                if rank([basis_SF', basis(j, :)']) ~= rank(basis_SF')
                
                    % Add this vector to basis_SF
                    basis_SF(end+1, :) = R(basis_reac_num(j), :);
                    
                    % Take note which reaction vectors from 'basis' are added to basis_SF
                    added_reaction(end+1) = j;
                    
                    % Stop forming basis_SF when its size reaches the rank s of the network
                    % Note: basis_SF will also have a rank s since it is composed of linearly independent vectors that span R
                    if size(basis_SF, 1) == size(basis, 1)
                        break
                    end
                end
            end
        end    
        
        % Form the second set of the separated basis vectors
        B2 = basis_reac_num(added_reaction);
        
        
        
    %
    % STEP 10: Check if R is the union of span(B1) and span(B2)
    %
        
        % Form span_B1 and span_B2
        % binary_decomp is 1 if R is the union of span(B1) and span(B2), 0 otherwise
        [binary_decomp, span_B1] = R_is_span_union(B1, B2, R);
        
        % Initialize list of power sets of B2 (for use in case we don't get an independent binary decomposition)
        power_set_B2 = { };
        
        % Create a list of power sets of B2
        % We exclude the empty set and the full set
        for j = 1:numel(B2)-1
            power_set_B2{j} = nchoosek(B2, j);
        end
        
        % Case 1: R is NOT the union of span(B1) and span(B2)
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
                    
                    % Check if R is the union of span(B1) and span(B2)
                    [binary_decomp, span_B1] = R_is_span_union(B1_new, B2_new, R);
                    
                    % Once successful in getting an independent binary decomposition
                    if binary_decomp == 1
                        
                        
                        
    %
    % STEP 12: Get the deficiency of the network induced by span(B1)
    %
                        
                        % Create network N1 called model_N1 from span(B1) and compute its deficiency delta1
                        [model_N1, delta1] = deficiency(model, span_B1);
                        
                        
                        
    %
    % STEP 13: Repeat from STEP 11 until the deficiency of the network induced by span(B1) is less than or equal to 1
    %
                        
                        % Case 1: The deficiency is 0
                        if delta1 == 0
                            
                            
                            
    %
    % STEP 14: For deficiency 0, check if the induced network is a weakly reversible power law system with reactant-determined kinetics (PL-RDK) with the Shinar-Feinberg pair in the same linkage class in the induced network
    %
                            
                            PL_RDK = is_PL_RDK(model_N1, m);
                            if PL_RDK
                                
                                % Check that model_N1 is weakly reversible
                                weakly_reversible = is_weakly_reversible(model_N1, m);
                                if weakly_reversible
                                    
                                    % Get the reactant complex of each reaction in the Shinar-Feinberg pair
                                    reactant_complex1 = reactant_complex_list{SF_pair1};
                                    reactant_complex2 = reactant_complex_list{SF_pair2};
                                    
                                    % Check if the reactant complexes are in the same linkage class in model_N1
                                    same_linkage_class = in_same_linkage_class(reactant_complex1, reactant_complex2, model_N1, m);
                                    
                                    % If the reactant complexes are in the same linkage class
                                    if same_linkage_class
                                        
                                        % Add the species with its associated Shinar-Feinberg pair to the list ACR_species
                                        ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                        
                                        % This is created so we can exit the iterated loops
                                        flag = 1;
                                        
                                        % This breaks out of the loop for k
                                        break
                                    end
                                end
                            end
                        
                        % Case 2: The deficiency is 1
                        elseif delta1 == 1
                             
                          
                          
    %
    % STEP 15: For deficiency 1, check if the induced network is a PL-RDK system AND that the reactant complexes of the Shinar-Feinberg pairs are nonterminal complexes
    %
    
                            PL_RDK = is_PL_RDK(model_N1, m);
                            if PL_RDK
                                
                                % Get the reactant complex of each reaction in the Shinar-Feinberg pair
                                reactant_complex1 = reactant_complex_list{SF_pair1};
                                reactant_complex2 = reactant_complex_list{SF_pair2};
                                
                                % Check that the reactant complexes are nonterminal
                                if (is_nonterminal(reactant_complex1, model_N1, m) && is_nonterminal(reactant_complex2, model_N1, m))
                                
                                    % Add the species with its associated Shinar-Feinberg pair to the list ACR_species
                                    ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                    
                                    % This is created so we can exit the iterated loops
                                    flag = 1;
                                    
                                    % This breaks out of the loop for k
                                    break
                                end
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
            [model_N1, delta1] = deficiency(model, span_B1);
            
            % STEP 13
            % Case 1: The deficiency is 0
            if delta1 == 0
                
                % STEP 14
                PL_RDK = is_PL_RDK(model_N1, m);
                if PL_RDK
                    weakly_reversible = is_weakly_reversible(model_N1, m);
                    if weakly_reversible
                        reactant_complex1 = reactant_complex_list{SF_pair1};
                        reactant_complex2 = reactant_complex_list{SF_pair2};
                        same_linkage_class = in_same_linkage_class(reactant_complex1, reactant_complex2, model_N1, m);
                        if same_linkage_class
                            ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                            % No break here since we don't want to exit checking the other pairs
                        end
                    end
                end

            % Case 2: The deficiency is 1
            elseif delta1 == 1
            
                % STEP 15
                PL_RDK = is_PL_RDK(model_N1, m);
                if PL_RDK
                    reactant_complex1 = reactant_complex_list{SF_pair1};
                    reactant_complex2 = reactant_complex_list{SF_pair2};
                    if (is_nonterminal(reactant_complex1, model_N1, m) && is_nonterminal(reactant_complex2, model_N1, m))
                        ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                        % No break here since we don't want to exit checking the other pairs
                    end
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
                        [binary_decomp, span_B1] = R_is_span_union(B1_new, B2_new, R);
                        if binary_decomp == 1
                            
                            % STEP 12
                            [model_N1, delta1] = deficiency(model, span_B1);
                            
                            % STEP 13
                            if delta1 == 0
                                
                                % STEP 14
                                PL_RDK = is_PL_RDK(model_N1, m);
                                if PL_RDK
                                    weakly_reversible = is_weakly_reversible(model_N1, m);
                                    if weakly_reversible
                                        reactant_complex1 = reactant_complex_list{SF_pair1};
                                        reactant_complex2 = reactant_complex_list{SF_pair2};
                                        same_linkage_class = in_same_linkage_class(reactant_complex1, reactant_complex2, model_N1, m);
                                        if same_linkage_class
                                            ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                            flag = 1;
                                            break
                                        end
                                    end
                                end
                            elseif delta1 == 1
                                
                                % STEP 15
                                PL_RDK = is_PL_RDK(model_N1, m);
                                if PL_RDK
                                    reactant_complex1 = reactant_complex_list{SF_pair1};
                                    reactant_complex2 = reactant_complex_list{SF_pair2};
                                    if (is_nonterminal(reactant_complex1, model_N1, m) && is_nonterminal(reactant_complex2, model_N1, m))
                                        ACR_species{end+1} = [model.species{SF_pair(i, 3)} ' [SF-pair (R' num2str(SF_pair(i, 1)) ', R' num2str(SF_pair(i, 2)) ') | deficiency ' num2str(delta1) ']'];
                                        flag = 1;
                                        break
                                    end
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

end










% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                                                   % %
% % The following are functions used in the algorithm % %
% %                                                   % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 1 of 15: model_species                                     %
%                                                                     %
%    - Purpose: To fill out list of species based on given reactions  %
%    - Input                                                          %
%         - model: empty species list                                 %
%    - Outputs                                                        %
%         - model: completed structure                                %
%         - m: number of species                                      %
%    - Used in                                                        %
%         - acr (STEP 1)                                              %
%         - deficiency                                                %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model, m] = model_species(model)

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
    
    % Count the number of species
    m = numel(model.species);
    
end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
% Function 2 of 15: stoich_matrix                                 %
%                                                                 %
%    - Purpose: To form the stoichometrix matrix N and the set of %
%          reactant and product complexes                         %
%    - Inputs                                                     %
%         - model: complete structure                             %
%         - m: number of species                                  %
%    - Outputs                                                    %
%         - N: stoichiometric matrix                              %
%         - reactant_complex: matrix of reactant complexes        %
%         - product_complex: matrix of product complexes          %
%         - r: total number of reactions                          %
%    - Used in                                                    %
%         - acr (STEP 2)                                          %
%         - deficiency                                            %
%         - is_PL_RDK                                             %
%         - is_weakly_reversible                                  %
%         - in_same_linkage_class                                 %
%         - is_nonterminal                                        %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [N, reactant_complex, product_complex, r] = stoich_matrix(model, m)

    % Initialize the matrix of reactant complexes
    reactant_complex = [ ];
    
    % Initialize the matrix of product complexes
    product_complex = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complex(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complex(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complex(:, end) - reactant_complex(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complex(:, end+1) = product_complex(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complex(:, end+1) = reactant_complex(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Count the total number of reactions
    r = size(N, 2);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                   %
% Function 3 of 15: kin_ord_matrix                                  %
%                                                                   %
%    - Purpose: To form the kinetic order matrix of the mass action %
%         system                                                    %
%    - Inputs                                                       %
%         - model: complete structure                               %
%         - m: number of species                                    %
%    - Output                                                       %
%         - F: kinetic order matrix                                 %
%    - Used in                                                      %
%         - acr (STEP 3)                                            %
%         - is_PL_RDK                                               %
%                                                                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function F = kin_ord_matrix(model, m)
    
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

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% Function 4 of 15: R_is_span_union                                     %
%                                                                       %
%    - Purpose: To check if R is the union of the span(B1) and span(B2) %
%    - Inputs                                                           %
%         - B1: set of vectors to be extended to a basis                %
%         - B2: set of added vectors to form a basis                    %
%         - R: reaction matrix                                          %
%    - Ouputs                                                           %
%         - binary_decomp: logical; whether R is the union of span(B1)  %
%              and span(B2) or not                                      %
%         - span_B1: span(B1)                                           %
%    - Used in acr (STEPS 10, 11)                                       %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [binary_decomp, span_B1] = R_is_span_union(B1, B2, R)
    
    % Get the reaction vectors forming B1 and forming B2
    set_B1 = R(B1, :);
    set_B2 = R(B2, :);
    
    % Initialize list of reactions in span(B1) and in span(B2)
    span_B1 = [ ];
    span_B2 = [ ];
    
    % Go through each reaction vector
    for i = 1:size(R, 1)
        
        % If the reaction vector is a linear combination of set_B1
        if rank([set_B1', R(i, :)']) == rank(set_B1')
            
            % Add this reaction vector number to span_B1
            span_B1(end+1) = i;
        end
        
        % If the reaction vector is a linear combination of set_B2
        if rank([set_B2', R(i, :)']) == rank(set_B2')
            
            % Add this reaction vector number to span_B2
            span_B2(end+1) = i;
        end
    end
    
    % Get the union of span(B1) and span(B2)
    span_B1_U_span_B2 = union(span_B1, span_B2);
    
    % If the union becomes a vertical vector: This happens in Octave when one of the sets being combined is empty
    if size(span_B1_U_span_B2, 1) > 1
        
        % Transpose
        span_B1_U_span_B2 = span_B1_U_span_B2';
    end
    
    % Check if span_B1_U_span_B2 contains all reactions
    % If R is the union of span(B1) and span(B2), an independent binary decomposition is formed
    if isequal(1:size(R, 1), span_B1_U_span_B2)
        binary_decomp = 1;
    else
        binary_decomp = 0;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 5 of 15: init_graph                                        %
%                                                                     %
%    - Purpose: To initialize an undirected graph                     %
%    - Input: none                                                    %
%    - Output                                                         %
%         - g: empty structure with subfields 'vertices' and 'edges'  %
%    - Used in                                                        %
%         - linkage_classes                                           %
%         - strong_linkage_classes                                    %
%         - in_same_linkage_class                                     %
%         - is_nonterminal                                            %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = init_graph()
    
    % Initialize a structure with empty fields
    g = struct();
    
    % Initialize 'vertices' field
    g.vertices = cell();
    
    % Initialize 'edges' field
    g.edges = cell();

end



% % % % % % % % % % % % % % % % % % % % % % % %
%                                             %
% Function 6 of 15: add_vertex                %
%                                             %
%    - Purpose: To add a vertex to a graph g  %
%    - Inputs                                 %
%         - g: graph structure                %
%         - v: name of vertex to be added     %
%    - Output                                 %
%         - g: structure with vertex added    %
%    - Used in                                %
%         - linkage_classes                   %
%         - strong_linkage_classes            %
%         - in_same_linkage_class             %
%         - is_nonterminal                    %
%                                             %
% % % % % % % % % % % % % % % % % % % % % % % %

function g = add_vertex(g, v)
    
    % Determine the index of v, if it already exists, in the 'vertices' field
    location = find(strcmp(v, g.vertices));
    
    % Case 1: 'location' is empty
    if isempty(location)
        
        % Add vertex v in the list of vertices in g
        g.vertices{end+1} = v;
        
        % Get the vertex number of the added vertex
        vertex_num = find(strcmp(v, g.vertices));
        
        % Initialize the place in g.edges where the edges connecting v to other vertices will be indicated
        g.edges{vertex_num} = [ ];
    
    % Case 2: 'location' is not empty
    else
        disp(['Vertex ' v ' is already in the graph.']);
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
% Function 7 of 15: add_edge                            %
%                                                       %
%    - Purpose: To add an undirected edge to a graph g  %
%    - Inputs                                           %
%         - g: graph structure with vertices v1 and v2  %
%         - v1: one vertex of edge to be added          %
%         - v2: another vertex of edge to be added      %
%    - Output                                           %
%         - g: structure with edge added                %
%    - Used in                                          %
%         - linkage_classes                             %
%         - in_same_linkage_class                       %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = add_edge(g, v1, v2)
        
    % Case 1: v1 or v2 is not in g
    if (isempty(find(strcmp(v1, g.vertices))) || isempty(find(strcmp(v2, g.vertices))))
        disp(['Make sure both ' v1 ' and ' v2 ' are in the graph.']);
        
    % Case 2: Both vertices are already in g
    else
        
        % Case 2.1: The vertices are the same
        if strcmp(v1, v2) == 1
            disp(['Make sure the vertices are different.']);
    
        % Case 2.2: The vertices are different
        else
            
            % Get the index of v1 and v2 in g.vertices
            v1_index = find(strcmp(v1, g.vertices));
            v2_index = find(strcmp(v2, g.vertices));
            
            % Check if an edge with v1 already exists
            try
                g.edges{v1_index};
            catch
                
                % If none, initialize g.edges for v1
                g.edges{v1_index} = cell();
            end
            
            % Check if an edge with v2 already exists
            try
                g.edges{v2_index};
            catch
                
                % If none, initialize g.edges for v2
                g.edges{v2_index} = cell();
            end
            
            % After all the controls above have been implemented, add an edge in the 'edges' field in g for both v1 and v2
            g.edges{v1_index}(end+1) = struct('vertex', v2_index, 'label', [v1 '-' v2]);
            g.edges{v2_index}(end+1) = struct('vertex', v1_index, 'label', [v2 '-' v1]);
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
% Function 8 of 15: linkage_classes                                         %
%                                                                           %
%    - Purpose: To determine the linkage class where each vertex belongs to %
%    - Inputs                                                               %
%         - r: number of reactions                                          %
%         - n: number of complexes                                          %
%         - Y: matrix of complexes                                          %
%         - reactant_complex: matrix of reactant complexes                  %
%         - product_complex: matrix of product complexes                    %
%         - model: complete structure                                       %
%    - Output                                                               %
%         - l: list of linkage classes where each vertex belongs to         %
%    - Used in                                                              %
%         - deficiency                                                      %
%         - is_weakly_reversible                                            %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function l = linkage_classes(r, n, Y, reactant_complex, product_complex, model)
    
    % Initialize an undirected graph g
    g = init_graph();

    % Go through each column of Y (a complex)
    for i = 1:n
        
        % For the zero complex
        if numel(find(Y(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(Y(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [model.species{find(Y(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(Y(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of g
        g = add_vertex(g, complex);
    end
    
    % Add edges to g: Ci -> Cj forms an edge
    % ~ suppresses the original output
    for i = 1:r
        g = add_edge(g, g.vertices{[~, loc] = ismember(reactant_complex(:, i)', Y', 'rows')}, g.vertices{[~, loc] = ismember(product_complex(:, i)', Y', 'rows')});
    end
    
    % Initialize the vector which will indicate in which linkage class number a vertex (i.e., complex) belongs to
    linkage_class = zeros(numel(g.vertices), 1);
    
    % Initialize the linkage class number tracker
    linkage_class_num = 0;
    
    % Go to each vertex
    for i = 1:numel(g.vertices)
        
        % Pay attention only to a vertex which has no linkage class number yet
        if linkage_class(i) == 0
            
            % This vertex goes to the next linkage class number
            linkage_class_num += 1;
            
            % Assign the linkage class number to the vertex
            linkage_class(i) = linkage_class_num;
            
            % Take note of the vertex that needs to be checked for edges
            to_check = [i];
            
            % Continue assigning a linkage class number to vertices that get into the check list
            while ~isempty(to_check)
                
                % Get the vertex in the check list
                v1 = to_check(end);
                
                % Remove the vertex from the check list (since we now know which vertex to focus on)
                to_check(end) = [ ];
                
                % Check to which vertices the vertex is connected to
                for j = 1:numel(g.edges{v1})
                    
                    % Take note of the vertex it is connected to
                    v2 = g.edges{v1}(j).vertex;
                    
                    % Pay attention to this vertex if it has no linkage class number yet
                    if linkage_class(v2) == 0
                        
                        % Assign this vertex with the same linkage class number as the vertex it is connected to
                        linkage_class(v2) = linkage_class_num;
                        
                        % Add this vertex to our check list: in the next round, we'll check to which other vertices it is connected to
                        to_check(end+1) = v2;
                    end
                end
            end
        end
    end
    
    % Count the number of linkage classes
    l = max(linkage_class);

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% Function 9 of 15: deficiency                                          %
%                                                                       %
%    - Purpose: To get the deficiency of a system created from span(B1) %
%    - Inputs                                                           %
%         - model: complete structure                                   %
%         - span_B1: span(B1)                                           %
%    - Ouputs                                                           %
%         - model_N1: new model created from span(B1)                   %
%         - delta1: deficiency of model_N1                              %
%    - Used in acr (STEP 12)                                            %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [model_N1, delta1] = deficiency(model, span_B1)
    
    % Name the new model
    model_N1.id = 'span_B1';
    
    % Prepare the list of species
    model_N1.species = { };
    
    % Create a vector of model.reaction numbers for the total number of reactions
    reac_num = [ ];
    for i = 1:numel(model.reaction)
        if model.reaction(i).reversible == 0
            reac_num(end+1) = i;
        else
            reac_num(end+1) = i;
            reac_num(end+1) = i;
        end
    end
    
    % Get all model reaction numbers corresponding to the reactions in span_B1
    reac = unique(reac_num(span_B1));
    
    % Put these reaction numbers in model_N1
    for i = 1:numel(reac)
        model_N1.reaction(i) = model.reaction(reac(i));
    end
    
    % Add to model_N1.species all species indicated in the reactions of model_N1
    [model_N1, m] = model_species(model_N1);
    
    % Form stoichiometric matrix N
    [N, reactant_complex, product_complex, r] = stoich_matrix(model_N1, m);
        
    % Get just the unique complexes
    % ind2(i) is the index in Y of the reactant complex in reaction i
    [Y, ~, ~] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
    
    % Count the number of linkage classes
    l = linkage_classes(r, n, Y, reactant_complex, product_complex, model);
    
    % Get the rank of the reaction network
    % S = Im N
    % dim S = dim (Im N) = rank(N)
    % Note: We talk of "dimension of a linear transformation" and "rank of a matrix"
    s = rank(N);
    
    % Compute the deficiency of the reaction network
    delta1 = n - l - s;

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
% Function 10 of 15: is_PL_RDK                                    %
%                                                                 %
%    - Purpose: To check if a system is a power law system with   %
%         reactant-determined kinetics (PL-RDK)                   %
%    - Inputs                                                     %
%         - model: complete structure                             %
%         - m: number of species                                  %
%    - Ouput                                                      %
%         - PL_RDK: logical; whether the system is PL-RDK or not  %
%    - Used in acr (STEPS 14, 15)                                 %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function PL_RDK = is_PL_RDK(model, m)
    
    % Get reactant_complex from function for stoichiometric matrix
    [~, reactant_complex, ~, ~] = stoich_matrix(model, m);
    
    % Form matrix of reactant complexes
    reactant_complex = reactant_complex';
    
    % Get the unique reactant complexes
    [reactant_complex_unique, ~, label] = unique(reactant_complex, 'rows');
    
    % Count the number of unique reactant complexes
    n_r = size(reactant_complex_unique, 1);
        
    % Get the labels of non-unique reactant complexes
    same_label = find(hist(label, unique(label)) > 1);
    
    % Initialize list reactions numbers with similar reactants
    branching_complex = { };
    
    % Group together reaction numbers with the same label
    for i = 1:numel(same_label)
        branching_complex{i} = find((label == same_label(i)));
    end
    
    % Initialize list of branching reactions
    branching_reaction = [ ];
    
    % Create a list of pairwise branching reactions
    for i = 1:numel(branching_complex)
        for j = 1:numel(branching_complex{i})
            for k = j+1:numel(branching_complex{i})
                
                % The list will have unique pairings: if (1, 2) is recorded, (2, 1) will no longer be noted
                branching_reaction(end+1, :) = [branching_complex{i}(j), branching_complex{i}(k)];
            end
        end
    end

    % Form kinetic order matrix F
    F = kin_ord_matrix(model, m);
    
    % Initialize list to track each pair
    F_identical = [ ];
    
    % Go through each pair
    for i = 1:size(branching_reaction, 1)
        
        % 1 means the pair has the same kinetic orders, 0 otherwise
        F_identical(i) = isequal(F(branching_reaction(i, 1), :), F(branching_reaction(i, 2), :));
    end
    
    % Case 1: There is at least 1 pair with different kinetic orders
    if ismember(0, F_identical)
        
        % The PLK system does NOT have reactant-determined kinetics (RDK)
        PL_RDK = 0;
    
    % Case 2: Each pair have the same kinetic order
    else
        
        % The PLK system has RDK
        PL_RDK = 1;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
% Function 11 of 15: add_path                           %
%                                                       %
%    - Purpose: To add a directed edge to a graph g     %
%    - Inputs                                           %
%         - g: graph structure with vertices v1 and v2  %
%         - v1: starting vertex of edge to be added     %
%         - v2: ending vertex of edge to be added       %
%    - Output                                           %
%         - g: structure with edge added                %
%    - Used in                                          %
%         - strong_linkage_classes                      %
%         - is_nonterminal                              %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = add_path(g, v1, v2)
    
    % Case 1: v1 or v2 is not in g
    if (isempty(find(strcmp(v1, g.vertices))) || isempty(find(strcmp(v2, g.vertices))))
        disp(['Make sure both ' v1 ' and ' v2 ' are in the graph.']);
        
    % Case 2: Both vertices are already in g
    else
        
        % Case 2.1: The vertices are the same
        if strcmp(v1, v2) == 1
            disp(['Make sure the vertices are different.']);
    
        % Case 2.2: The vertices are different
        else
            
            % Get the index of v1 and v2 in g.vertices
            v1_index = find(strcmp(v1, g.vertices));
            v2_index = find(strcmp(v2, g.vertices));
            
            % Check if an edge from v1 to v2 already exists
            for i = 1:numel(g.edges{v1_index})
                if g.edges{v1_index}(i).vertex == v2_index
                    disp(['Edge ' v1 '->' v2 ' already exists.']);
                    
                    % 'return' exits the function; we don't need to continue the code
                    % If we wanted to just get out of the loop, we use 'break'
                    return
                end
            end
            
            % Add a directed edge in the 'edges' field in g for v1
            g.edges{v1_index}(end+1) = struct('vertex', v2_index, 'label', [v1 '->' v2]);
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                                   %
% Function 12 of 15: strong_linkage_classes                                         %
%                                                                                   %
%    - Purpose: To determine the strong linkage class where each vertex belongs to  %
%    - Inputs                                                                       %
%         - r: number of reactions                                                  %
%         - n: number of complexes                                                  %
%         - Y: matrix of complexes                                                  %
%         - reactant_complex: matrix of reactant complexes                          %
%         - product_complex: matrix of product complexes                            %
%         - model: complete structure                                               %
%    - Output                                                                       %
%         - sl: list of strong linkage classes where each vertex belongs to         %
%    - Used in                                                                      %
%         - is_weakly_reversible                                                    %
%                                                                                   %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function sl = strong_linkage_classes(r, n, Y, reactant_complex, product_complex, model)
    
    % Initialize a directed graph g
    g = init_graph();

    % Go through each column of Y (a complex)
    for i = 1:n
        
        % For the zero complex
        if numel(find(Y(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(Y(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [model.species{find(Y(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(Y(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of G
        g = add_vertex(g, complex);
    end
    
    % Add a directed edge to g: Ci -> Cj forms an edge
    for i = 1:r
        g = add_path(g, g.vertices{[~, loc] = ismember(reactant_complex(:, i)', Y', 'rows')}, g.vertices{[~, loc] = ismember(product_complex(:, i)', Y', 'rows')});
    end
    
    % Define function which visits each complex (i.e., vertex) v and the other vertices connected to it
    function visit(v)
        
        % Set the discovery time of the complex as the current time
        discovery_time(v) = time;
        
        % Set the discovery time of the strong linkage class as the current time
        slc_discovery_time(v) = time;
        
        % Move the time forward for the next complex
        time = time + 1;
        
        % Add the complex in the list of complexes in the same strong linkage class
        stack(end+1) = v;
        
        % Note that the complex is already listed in the strong linkage class
        on_stack(v) = true;
        
        % Go through each edge connected to the vertex (i.e., complex)
        for j = 1:numel(g.edges{v})
            
            % Take the vertex connected to it
            v2 = g.edges{v}(j).vertex;
            
            % If the vertex is not yet visited
            if discovery_time(v2) == 0
                
                % Apply the visit function to this vertex
                visit(v2);
                
                % Set the discovery time of the strong linkage class
                % slc_discovery_time(v2) < slc_discovery_time(v) iff a vertex in the stack before v is reachable from v2 (and so they are all in the same strong linkage class)
                slc_discovery_time(v) = min(slc_discovery_time(v), slc_discovery_time(v2));
            
            % If v2 was visited before v
            elseif on_stack(v2)
                
                % So v and v2 are in the same component, and they must have the same slc_discovery_time
                slc_discovery_time(v) = min(slc_discovery_time(v), slc_discovery_time(v2));
            end
        end
        
        % If v is the first visited node of its strong linkage class, all the other vertices of the strong linkage class follow it on the stack
        if slc_discovery_time(v) == discovery_time(v) 
            while true
                v2 = stack(end);
                stack(end) = [ ];
                on_stack(v2) = false;
                strong_linkage_class(v2) = slc_num;
                if v2 == v
                    break
                end
            end
            
            % Add 1 to the strong linkage class number
            slc_num = slc_num + 1;
        end
    end
    
    % This is the actual function
    
    % Initialize the discovery time of the vertices
    discovery_time = zeros(numel(g.vertices), 1);
    
    % Initialize the discovery time of the [first visited complex of the] strong linkage class of the vertices
    slc_discovery_time = zeros(numel(g.vertices), 1);
    
    % Start the timer at 1
    time = 1;
    
    % Initialize the list of complexes in the strong linkage class
    stack = [ ];
    
    % Initialize that no vertex is on a strong linkage class
    on_stack = false(numel(g.vertices), 1);
    
    % Initialize vector of strong linkage class numbers
    strong_linkage_class = zeros(numel(g.vertices), 1);
    
    % Start numbering the strong linkage class at 1
    slc_num = 1;
    
    % Go through each complex (i.e., vertex)
    for i = 1:numel(g.vertices)
        
        % Use the visit function if it has no strong linkage classs number
        if strong_linkage_class(i) == 0
            visit(i);
        end
    end
    
    % Count the number of strong linkage classes
    sl = max(strong_linkage_class);

end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 13 of 15: is_weakly_reversible                             %
%                                                                     %
%    - Purpose: To check if a system is weakly reversible             %
%    - Inputs                                                         %
%         - model: complete structure                                 %
%         - m: number of species                                      %
%    - Ouput                                                          %
%         - weakly_reversible: logical; whether the system is weakly  %
%              reversible or not                                      %
%    - Used in acr (STEP 14)                                          %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function weakly_reversible = is_weakly_reversible(model, m)
    
    % Form stoichiometric matrix N
    [N, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    % Get just the unique complexes
    % ind2(i) is the index in Y of the reactant complex in reaction i
    [Y, ~, ind2] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
    
    % Count the number of linkage classes
    l = linkage_classes(r, n, Y, reactant_complex, product_complex, model);
    
    % Count the number of strong linkage classes
    sl = strong_linkage_classes(r, n, Y, reactant_complex, product_complex, model);

    % Weakly reversible if the number of linkage classes and the number of strong linkage classes are the same
    if sl == l
        weakly_reversible = 1;
    else
        weakly_reversible = 0;
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Function 14 of 15: in_same_linkage_class                                %
%                                                                         %
%    - Purpose: To check if two complexes are in the same linkage class   %
%    - Inputs                                                             %
%         - complex1: first complex                                       %
%         - complex2: second complex                                      %
%         - model: complete structure                                     %
%         - m: number of species                                          %
%    - Ouput                                                              %
%         - same_linkage_class: logical; whether the two complexes are in %
%              the same linkage class or not                              %
%    - Used in acr (STEP 14)                                              %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function same_linkage_class = in_same_linkage_class(complex1, complex2, model, m)
    
    % Form stoichiometric matrix N
    [N, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    % Get just the unique complexes
    [Y, ~, ~] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
        
    % Initialize an undirected graph g
    g = init_graph();

    % Go through each column of Y (a complex)
    for i = 1:n
        
        % For the zero complex
        if numel(find(Y(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(Y(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [model.species{find(Y(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(Y(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of g
        g = add_vertex(g, complex);
    end
    
    % Check if complex1 is a valid string of complexes
    if isempty(find(strcmp(g.vertices, complex1)))
        disp([complex1 ' is not a complex in the network.']);
    end
    
    % Check if complex2 is a valid string of complexes
    if isempty(find(strcmp(g.vertices, complex2)))
        disp([complex2 ' is not a complex in the network.']);
    end
    
    % If complex1 or complex2 is not valid, exit the function
    if (isempty(find(strcmp(g.vertices, complex1))) || isempty(find(strcmp(g.vertices, complex2))))
        same_linkage_class = [ ];
        return
    end
    
    % Add edges to g: Ci -> Cj forms an edge
    % ~ suppresses the original output
    for i = 1:r
        g = add_edge(g, g.vertices{[~, loc] = ismember(reactant_complex(:, i)', Y', 'rows')}, g.vertices{[~, loc] = ismember(product_complex(:, i)', Y', 'rows')});
    end
    
    % Initialize the vector which will indicate in which linkage class number a vertex (i.e., complex) belongs to
    linkage_class = zeros(numel(g.vertices), 1);
    
    % Initialize the linkage class number tracker
    linkage_class_num = 0;
    
    % Go to each vertex
    for i = 1:numel(g.vertices)
        
        % Pay attention only to a vertex which has no linkage class number yet
        if linkage_class(i) == 0
            
            % This vertex goes to the next linkage class number
            linkage_class_num += 1;
            
            % Assign the linkage class number to the vertex
            linkage_class(i) = linkage_class_num;
            
            % Take note of the vertex that needs to be checked for edges
            to_check = [i];
            
            % Continue assigning a linkage class number to vertices that get into the check list
            while ~isempty(to_check)
                
                % Get the vertex in the check list
                v1 = to_check(end);
                
                % Remove the vertex from the check list (since we now know which vertex to focus on)
                to_check(end) = [ ];
                
                % Check to which vertices the vertex is connected to
                for j = 1:numel(g.edges{v1})
                    
                    % Take note of the vertex it is connected to
                    v2 = g.edges{v1}(j).vertex;
                    
                    % Pay attention to this vertex if it has no linkage class number yet
                    if linkage_class(v2) == 0
                        
                        % Assign this vertex with the same linkage class number as the vertex it is connected to
                        linkage_class(v2) = linkage_class_num;
                        
                        % Add this vertex to our check list: in the next round, we'll check to which other vertices it is connected to
                        to_check(end+1) = v2;
                    end
                end
            end
        end
    end
    
    % 1 if they are in the same linkage class, 0 otherwise
    same_linkage_class = linkage_class(find(strcmp(g.vertices, complex1))) == linkage_class(find(strcmp(g.vertices, complex2)));

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                       %
% Function 15 of 15: is_nonterminal                                     %
%                                                                       %
%    - Purpose: To check if the complexes are nonterminal               %
%    - Inputs                                                           %
%         - check_complex: to check if nonterminal                      %
%         - model: complete structure                                   %
%         - m: number of species                                        %
%    - Ouput                                                            %
%         - nonterminal: logical; whether the complex is nonterminal or %
%              not                                                      %
%    - Used in acr (STEP 15)                                            %
%                                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function nonterminal = is_nonterminal(check_complex, model, m)

    % Form stoichiometric matrix N
    [N, reactant_complex, product_complex, r] = stoich_matrix(model, m);
    
    % Get just the unique complexes
    [Y, ~, ~] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
        
    % Initialize a directed graph g
    g = init_graph();

    % Go through each column of Y (a complex)
    for i = 1:n
        
        % For the zero complex
        if numel(find(Y(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(Y(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [model.species{find(Y(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(Y(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of G
        g = add_vertex(g, complex);
    end
    
    % Check if check_complex is a valid string of complexes
    if isempty(find(strcmp(g.vertices, check_complex)))
        disp([check_complex ' is not a complex in the network.']);
        nonterminal = [ ];
        return
    end
    
    % Add a directed edge to g: Ci -> Cj forms an edge
    for i = 1:r
        g = add_path(g, g.vertices{[~, loc] = ismember(reactant_complex(:, i)', Y', 'rows')}, g.vertices{[~, loc] = ismember(product_complex(:, i)', Y', 'rows')});
    end
    
    % Define function which visits each complex (i.e., vertex) v and the other vertices connected to it
    function visit(v)
        
        % Set the discovery time of the complex as the current time
        discovery_time(v) = time;
        
        % Set the discovery time of the strong linkage class as the current time
        slc_discovery_time(v) = time;
        
        % Move the time forward for the next complex
        time = time + 1;
        
        % Add the complex in the list of complexes in the same strong linkage class
        stack(end+1) = v;
        
        % Note that the complex is already listed in the strong linkage class
        on_stack(v) = true;
        
        % Go through each edge connected to the vertex (i.e., complex)
        for j = 1:numel(g.edges{v})
            
            % Take the vertex connected to it
            v2 = g.edges{v}(j).vertex;
            
            % If the vertex is not yet visited
            if discovery_time(v2) == 0
                
                % Apply the visit function to this vertex
                visit(v2);
                
                % Set the discovery time of the strong linkage class
                % slc_discovery_time(v2) < slc_discovery_time(v) iff a vertex in the stack before v is reachable from v2 (and so they are all in the same strong linkage class)
                slc_discovery_time(v) = min(slc_discovery_time(v), slc_discovery_time(v2));
            
            % If v2 was visited before v
            elseif on_stack(v2)
                
                % So v and v2 are in the same component, and they must have the same slc_discovery_time
                slc_discovery_time(v) = min(slc_discovery_time(v), slc_discovery_time(v2));
            end
        end
        
        % If v is the first visited node of its strong linkage class, all the other vertices of the strong linkage class follow it on the stack
        if slc_discovery_time(v) == discovery_time(v) 
            while true
                v2 = stack(end);
                stack(end) = [];
                on_stack(v2) = false;
                strong_linkage_class(v2) = slc_num;
                if v2 == v
                    break
                end
            end
            
            % Add 1 to the strong linkage class number
            slc_num = slc_num + 1;
        end
    end
    
    % This is the actual function
    
    % Initialize the discovery time of the vertices
    discovery_time = zeros(numel(g.vertices), 1);
    
    % Initialize the discovery time of the [first visited complex of the] strong linkage class of the vertices
    slc_discovery_time = zeros(numel(g.vertices), 1);
    
    % Start the timer at 1
    time = 1;
    
    % Initialize the list of complexes in the strong linkage class
    stack = [ ];
    
    % Initialize that no vertex is on a strong linkage class
    on_stack = false(numel(g.vertices), 1);
    
    % Initialize vector of strong linkage class numbers
    strong_linkage_class = zeros(numel(g.vertices), 1);
    
    % Start numbering the strong linkage class at 1
    slc_num = 1;
    
    % Go through each complex (i.e., vertex)
    for i = 1:numel(g.vertices)
        
        % Use the visit function if it has no strong linkage classs number
        if strong_linkage_class(i) == 0
            visit(i);
        end
    end
    
    % Initialize list of non-terminal strong linkage classes
    nonterminal_complexes = [ ];
    
    % Go through each complex
    for i = 1:numel(g.edges)
        
        % Check to which complex it is connected to
        for j = 1:numel(g.edges{i})
            
            % Check if the two complexes belong to the same strong linkage class
            same_strong_linkage_class = strong_linkage_class(i) == strong_linkage_class([g.edges{i}.vertex](j));
            
            % If they do not, then the complex does not belong in a terminal strong linkage class
            if same_strong_linkage_class == 0
                nonterminal_complexes(end+1) = i;
            end
        end
    end
    
    % Generate the list of complexes that do not belong to a terminal strong linkage class
    nonterminal_complexes = unique(nonterminal_complexes);
    
    % Locate the other complexes belonging to the same strong linkage class as the ones identified to not belong in a terminal strong linkage class
    % These constitute the nonterminal complexes
    nonterminal_complexes = find(ismember(strong_linkage_class, strong_linkage_class(nonterminal_complexes)));
    
    % Check if the input complex is nonterminal
    if isempty(find(strcmp(g.vertices(nonterminal_complexes), check_complex))) == 1
        nonterminal = 0;
    else
        nonterminal = 1;
    end

end