% Read the Excel file
fileName = './GMSMS.xlsx';
inputData = readtable(fileName);

% Select columns 2-4 from inputData
inputData = inputData{:, 2:4};

% 1. Calculate mean and variance for each row of inputData
means = mean(inputData, 2);
variances = var(inputData, 0, 2);

% 2. Generate random values based on mean and variance for each row
% 3. Save generated random values in Excel spreadsheet
% 4. Generate 3 groups of random values and save to different Excel files

% Set parameters
num_groups = 2;
num_random_values = 4;  % Number of random values per group, can be modified as needed

% Generate 3 groups of random values and save to different Excel files
for group = 1:num_groups
    random_values = zeros(size(inputData, 1), num_random_values);
    for i = 1:size(inputData, 1)
        random_values(i, :) = normrnd(means(i), sqrt(variances(i)), 1, num_random_values);
    end
    
    % Create table and set column names
    T = array2table(random_values, 'VariableNames', arrayfun(@(x) sprintf('Random_%d', x), 1:num_random_values, 'UniformOutput', false));
    
    % Save current group results to separate Excel file
    output_filename = sprintf('Random_Group%d_%dvalues_output.xlsx', group, num_random_values);
    writetable(T, output_filename);
    fprintf('Group %d results saved to %s\n', group, output_filename);
end