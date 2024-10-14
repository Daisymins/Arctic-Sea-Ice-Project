%% LOAD CSV DATA
file_names = [ ...
    "N_09_extent_v3.0.csv", ...
    "SIA_CMIP6_historical_September.txt", ...
    "SIA_CMIP6_ssp126_September.txt", ...
    "SIA_CMIP6_ssp245_September.txt", ...
    "SIA_CMIP6_ssp585_September.txt" ...
];

% Read all tables in a loop
tables = cell(1, length(file_names));
for i = 1:length(file_names)
    tables{i} = readtable(file_names{i}, 'VariableNamingRule', 'preserve');
end

% Assign tables to variables for easier reference
obsv_Septdata = tables{1};
historical_Sept = tables{2};
ssp126_model = tables{3};
ssp245_model = tables{4};
ssp585_model = tables{5};

%% ALIGN DATA (1979-2014 common years)
obs_years = obsv_Septdata.year;
obs_area = obsv_Septdata.area;
model_years_hist = historical_Sept.Year;

% Find common years and extract relevant data
common_years = intersect(obs_years, model_years_hist);
obs_common = obs_area(ismember(obs_years, common_years));
model_common_hist = historical_Sept(ismember(model_years_hist, common_years), :);

%% Calculate Observed Data Mean and Confidence Interval
obs_mean = mean(obs_common);
obs_std = std(obs_common);
n_obs = length(obs_common);
z = 1.96; % 95% confidence level

% Confidence interval calculation
ci_lower = obs_mean - z * (obs_std / sqrt(n_obs));
ci_upper = obs_mean + z * (obs_std / sqrt(n_obs));

disp(['Observed Mean: ', num2str(obs_mean)]);
disp(['95% Confidence Interval: [', num2str(ci_lower), ', ', num2str(ci_upper), ']']);

%% Model Selection: Calculate Model Means and Compare with Confidence Interval
model_names_hist = historical_Sept.Properties.VariableNames(2:end);
n_models_hist = length(model_names_hist);
model_means = zeros(1, n_models_hist);
is_within_ci = false(1, n_models_hist);

for i = 1:n_models_hist
    model_data_hist = model_common_hist{:, i + 1}; % +1 because first column is 'Year'
    model_mean = mean(model_data_hist);
    model_means(i) = model_mean;
    is_within_ci(i) = model_mean >= ci_lower && model_mean <= ci_upper; % Check within CI
end

% Display model comparison results
model_ci_comparison = table(model_names_hist', model_means', is_within_ci', ...
    'VariableNames', {'Model', 'Model_Mean', 'Within_CI'});
disp('Model Mean Comparison with Observed Confidence Interval:');
disp(model_ci_comparison);

%% Step 1: Perform Error Analysis on Historical Data (for Selected Models)
selected_models = {'NESM3', 'MRI-ESM2-0', 'MIROC-ES2L', 'CNRM-ESM2-1'};
valid_models_data_hist = historical_Sept(:, selected_models);

% Reuse common_years and observed data
obs_common = obs_area(ismember(obs_years, common_years));
model_common_hist = valid_models_data_hist(ismember(model_years_hist, common_years), :);

% Initialize error metric arrays
n_models_hist = numel(selected_models);
bias_hist = zeros(1, n_models_hist);
var_hist = zeros(1, n_models_hist);
std_dev_hist = zeros(1, n_models_hist);

% Consolidate error analysis into a single loop
for i = 1:n_models_hist
    model_data_hist = table2array(model_common_hist(:, i));
    bias_hist(i) = mean(model_data_hist - obs_common);
    var_hist(i) = var(model_data_hist);
    std_dev_hist(i) = std(model_data_hist);
    
    % Display error metrics
    fprintf('Error Analysis for Model %s:\n', selected_models{i});
    fprintf('  Bias: %.3f\n', bias_hist(i));
    fprintf('  Variance: %.3f\n', var_hist(i));
    fprintf('  Standard Deviation: %.3f\n\n', std_dev_hist(i));
end

%% Step 2: Normalize the Error Metrics
norm_var = (var_hist - min(var_hist)) / (max(var_hist) - min(var_hist));
norm_sd = (std_dev_hist - min(std_dev_hist)) / (max(std_dev_hist) - min(std_dev_hist));
norm_bias = (bias_hist - min(bias_hist)) / (max(bias_hist) - min(bias_hist));

%% Step 3: Define Weights Based on Error Metrics
w_bias = 0.5; w_var = 0.3; w_sd = 0.2;    
weights = 1 ./ (w_bias * norm_bias + w_var * norm_var + w_sd * norm_sd);
normalized_weights = weights / sum(weights);

% Debugging: Display normalized weights
%disp('Normalized Weights (sum to 1):');
%disp(normalized_weights);

%% Combine Model Predictions for Different SSP Scenarios



% Combine for each SSP scenario
[sea_ice_126_combined, years_126] = combine_model_predictions(ssp126_model, selected_models, normalized_weights);
[sea_ice_245_combined, years_245] = combine_model_predictions(ssp245_model, selected_models, normalized_weights);
[sea_ice_585_combined, years_585] = combine_model_predictions(ssp585_model, selected_models, normalized_weights);

%% Plot Combined Sea Ice Area Predictions for Different SSP Scenarios
figure;
plot(years_126, sea_ice_126_combined, '-b', 'LineWidth', 2); hold on;
plot(years_245, sea_ice_245_combined, '-r', 'LineWidth', 2);
plot(years_585, sea_ice_585_combined, '-g', 'LineWidth', 2);
xlabel('Year');
ylabel('Sea Ice Area (million km²)');
title('Combined Sea Ice Area Predictions for Different SSP Scenarios');
legend('SSP126 Combined', 'SSP245 Combined', 'SSP585 Combined');
grid on;

%% Find ice-free years for each SSP scenario
ice_free_threshold = 1;


% Find ice-free year for combined SSP126 scenario
ice_free_year_126_combined = find_ice_free_year(years_126, sea_ice_126_combined, ice_free_threshold);
if ~isempty(ice_free_year_126_combined)
    fprintf('The Arctic becomes ice-free under the SSP126 scenario in the year %d (combined models).\n', (ice_free_year_126_combined));
else
    disp('The Arctic does not become ice-free under the SSP126 scenario in the available data (combined models).');
end

% Find ice-free year for combined SSP245 scenario
ice_free_year_245_combined = find_ice_free_year(years_245, sea_ice_245_combined, ice_free_threshold);
if ~isempty(ice_free_year_245_combined)
    fprintf('The Arctic becomes ice-free under the SSP245 scenario in the year %d (combined models).\n', (ice_free_year_245_combined));
else
    disp('The Arctic does not become ice-free under the SSP245 scenario in the available data (combined models).');
end

% Find ice-free year for combined SSP585 scenario
ice_free_year_585_combined = find_ice_free_year(years_585, sea_ice_585_combined, ice_free_threshold);
if ~isempty(ice_free_year_585_combined)
    fprintf('The Arctic becomes ice-free under the SSP585 scenario in the year %d (combined models).\n', (ice_free_year_585_combined));
else
    disp('The Arctic does not become ice-free under the SSP585 scenario in the available data (combined models).');
end

%% Visualize Historical Error Analysis Results
figure;
bar(categorical(selected_models), [bias_hist', var_hist', std_dev_hist'], 'grouped');
xlabel('Model');
ylabel('Error Metric Value');
title('Error Analysis (Bias, Variance, Standard Deviation) for Selected Models (Historical Data)');
legend('Bias', 'Variance', 'Standard Deviation');
grid on;

%% Functions

% Function to find the ice-free year
function ice_free_year = find_ice_free_year(years, sea_ice_data, threshold)
    row = find(sea_ice_data < threshold, 1); % Find the first row where sea ice is below the threshold
    if ~isempty(row)
        ice_free_year = years(row); % Return the corresponding year
    else
        ice_free_year = []; % Return an empty array if no ice-free year is found
    end
end


% Create a function to combine model predictions based on weights
function [sea_ice_combined, years] = combine_model_predictions(model, selected_models, normalized_weights)
    sea_ice_combined = zeros(height(model), 1);
    valid_weights_sum = 0;  % Sum of weights for valid models
    
    for i = 1:numel(selected_models)
        model_data = model.(selected_models{i});  % Extract model data
        
        % Check if the model data contains NaN values
        if all(~isnan(model_data))  % Only include if there are no NaNs
            sea_ice_combined = sea_ice_combined + normalized_weights(i) * model_data;
            valid_weights_sum = valid_weights_sum + normalized_weights(i);  % Sum valid weights
        end
    end
    
    % Normalize the combined output by the sum of valid weights
    if valid_weights_sum > 0  % Avoid division by zero
        sea_ice_combined = sea_ice_combined / valid_weights_sum;
    else
        warning('No valid models found for combining predictions.');
    end
    
    years = model.Year; % Return the years associated with the combined sea ice area
end