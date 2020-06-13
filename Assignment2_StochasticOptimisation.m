[costs, probabilityVec] = find_best(1000, 500, 0.05);

%Plots the costs
plot(costs);

%Print the probability vector clamped to the range 0-1
min(max(probabilityVec, 0), 1)

%Gets the total cost for the last iteration
function [costs, probabilityVec] = find_best(iterations, solutionsPerEpoch, learnRate)
    probabilityVec = zeros(12, 3) + 0.5; %Creates a probability vector value with equal probability
    costs = linspace(0, 0, iterations); %Creates an array of cost values per each iteration
    lastCost = 9999;
    for i = 1:iterations
        [binaryVec, thisCost] = run_epoch(probabilityVec, solutionsPerEpoch);

        %If the current cost is better than the last cost, mutate
        if lastCost > thisCost  
            mutations = (binaryVec + (binaryVec - 1)) * learnRate; 

            %Mutate the probability vector minorly to increase likelihood of good values increasing
            probabilityVec = probabilityVec + mutations;
            lastCost = thisCost;
        end
        costs(i) = lastCost; %Save the cost at the iteration
    end
end

%Find the optimal binary vector for an input probability vector
function [binaryVec, cost] = run_epoch(probabilityVec, num_samples)

    binaryVec = probabilityVec < rand(size(probabilityVec));
    cost = calculate_cost(binaryVec);

    %Loop through the number of samples
    for i = 1:num_samples - 1
        newVec = probabilityVec > rand(size(probabilityVec)); 
        newCost = calculate_cost(newVec);
        if newCost < cost %Check if the newCost is less than the current cost. If so, save the cost
            cost = newCost;
            binaryVec = newVec;
        end
    end
end

%Calculate the cost of the epoch of sample solutions using cost table - from the inputted probability vector
function [cost] = calculate_cost(binaryVec)

    cost = 0; %Initialise cost to be 0

    %Given cost table (editable) - 8x12 matrix
    costTable = [15, 24, 30, 31, 26, 32, 23, 37, 32, 20, 13, 16;
                27, 18, 19, 36, 12, 17, 40, 33, 18, 39, 26, 32;
                37, 13, 15, 26, 19, 17, 36, 35, 30, 10, 25, 23;
                34, 25, 24, 32, 16, 35, 10, 38, 26, 39, 13, 14;
                32, 25, 15, 10, 36, 19, 33, 16, 11, 32, 15, 12;
                34, 18, 11, 30, 39, 35, 37, 15, 38, 18, 18, 15;
                20, 23, 11, 39, 34, 28, 36, 28, 32, 28, 25, 32;
                17, 34, 35, 40, 35, 28, 25, 28, 12, 13, 40, 22];

    terminals_per_concentrator = linspace(0, 0, 8); %Initialises how many terminals are connected to each concentrator

    
    for terminal = 1:size(binaryVec, 1)
        concentrator = binaryToDecimal(binaryVec(terminal, :)) + 1; %Convert the concentrator number to a decimal value
        terminals_per_concentrator(concentrator) = terminals_per_concentrator(concentrator) + 1; 
        cost = cost + costTable(concentrator, terminal);
    end
    
    %If there's any more than 3 terminals per concentrator, disregard this solution
    if max(terminals_per_concentrator) > 3
        cost = 999999;
    end
end

%Converts binary matrix to decimal
 function decimal = binaryToDecimal(binaryVec)
    decimal = binaryVec(3)*4 + binaryVec(2)*2 + binaryVec(1);
 end