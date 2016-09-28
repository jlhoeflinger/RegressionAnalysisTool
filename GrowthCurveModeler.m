function GrowthCurveModeler( file_or_dir, varargin)

%GrowthCurveModeler - Calculates some metrics for growth curves, as well as graphing
%a regression curve of the data points. 
%
%Required parameter: 
% file_or_dir - input data file or directory of input files, must be formatted properly
%
% Optional parameters:
%
% MaxTimepoint - value {default is total specified in dataset}
%   Final timepoint of the dataset 
%
% Threshold - value {default 0.3}
%   Threshold which describes the minimum OD reading to signify growth (default 0.3)
% 
% Model - [{'modlogistic'} | 'gompertz' | 'logistic' | 'modgompertz']
%   Regression model to plot the curves
%
% IncubationTime - value {default 1.0}
%   If you incubate your cells before recording timepoints,
%   this will appropriately shift your data points in time
%   for lag time calculations 
%
% DoubleHump - ['on' | {'off'}] 
%   Flag to indicate double hump processing, this expects all datasets to include a
%   double/multi hump and should remove all growth humps before the "main curve" 
% 
%
%
%  Examples:
%     GrowthCurveModeler('dataset.xlsx', 'DoubleHump', 'on', 'Threshold', 0.2);
%
%     GrowthCurveModeler('.', 'IncubationTime', 1.5);
% 
%     GrowthCurveModeler('folder_containing_xlsxfiles');




if (isdir(file_or_dir))
    runnable_files = dir([file_or_dir '/' '*.xlsx']);

    for i = 1:length(runnable_files)
        GrowthCurveModeler([file_or_dir '/' runnable_files(i).name], varargin{:}); 
    end
    return;
end

max_timepoint = -1;
growth_threshold = 0.3;
model = 'modlogistic';
incubation_time = 1.0;
double_hump = 0; 

i = 1;
while (i < length(varargin))    
   if (strcmp(varargin{i}, 'MaxTimepoint'))
       i = i + 1;
       max_timepoint = varargin{i};
   elseif (strcmp(varargin{i}, 'Threshold'))
       i = i + 1;
       growth_threshold = varargin{i};
   elseif (strcmp(varargin{i}, 'IncubationTime'))
       i = i + 1;
       incubation_time = varargin{i};
   elseif (strcmp(varargin{i}, 'DoubleHump'))
       i = i + 1;
       if(strcmp(varargin{i},'on'))
           double_hump = 1;
       end % else off
   elseif (strcmp(varargin{i}, 'Model'))
       i = i + 1;
       model = varargin{i};
   end  
   i = i+1;
end
    


[path, filestub, ext] = fileparts(file_or_dir); 

[Data, title_data] = xlsread(file_or_dir);
dims = size(title_data);



if (~isempty(path))
    path = [path '/'];
end
plots_folder = [path 'results/' filestub ' plots/'];
        
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end

output(1,:) = {'Sugar', 'Strain', 'Lag Time (hours)', 'Max Specific Growth Rate (1/hours)',  'Doubling Time (hours)', 'Max OD', 'Median OD', 'Delta OD', 'Notes', 'R^2', 'SSE' , 'RMSE'};
time_interval =0.5;

sugar = '';
strain_count = 0;
sugar_count = 0;

first_sugar = 1;
Sugars = {};
Start_idxs = [];
Strain_counts = [];
lag_times = [];
for i=1:dims(2)
    if (~isempty(char(title_data(1,i))))
        if  (~first_sugar)
            Strain_counts(sugar_count) = strain_count;
        end
        strain_count = 0;
        sugar_start_idx = i;
        first_sugar = 0;
        sugar_count = sugar_count + 1;
        
        sugar = char(title_data(1,i));
        Start_idxs(sugar_count) = sugar_start_idx;
        Sugars(sugar_count) = {char(sugar)};

    end
    if (strcmpi(char(title_data(2,i)),'Time'))
        time_interval = Data(2,i) - Data(1,i);
    else
        strain_count = strain_count + 1;
        strain = char(title_data(2,i));
        h = figure();set(h, 'WindowStyle', 'docked');
        if (max_timepoint < 0)
            [lagtime, max_u, OD_max, median_OD_max, delta_OD_max, doubling_time, note, goodness] = MicrobialKinetics(Data(:,i), time_interval, incubation_time, growth_threshold, model, double_hump);
        else
            [lagtime, max_u, OD_max, median_OD_max, delta_OD_max, doubling_time, note, goodness] = MicrobialKinetics(Data(1:max_timepoint/time_interval,i), time_interval, incubation_time, growth_threshold, model, double_hump);
        end
        lag_times(i) = lagtime;
        
        output(i,:) = {sugar,strain, lagtime, max_u, doubling_time, OD_max, median_OD_max, delta_OD_max, note, goodness.rsquare,  goodness.sse, goodness.rmse};
        name = [sugar '-' strain];
        title( name );
        sugar_folder = [plots_folder '/' sugar];
        if ~exist(sugar_folder, 'dir')
            mkdir(sugar_folder);
        end
        saveas(h,[sugar_folder '/' name], 'bmp');
        close(h);
    end
    
end

output_file = [ path 'results/' filestub ' results.xlsx'];
if (exist(output_file, 'file'))
	delete(output_file);
end
xlswrite(output_file, output);



end



function [ lag_time, max_spec_growth_rate, max_od,  median_od_max, delta_OD_max, doubling_time, note, goodness ] = MicrobialKinetics(OD_values, time_interval, incubation_time, threshold, model, double_hump)
%MicrobialKinetics -  For a specific dataset, attempts to determine a group
%of statistics on the growth which occurred.   Plots the dataset along the 
%way, with the best fit regression.



%  time interval assumed to be half hour time blocks

%set to zero initially

[max_od, max_location] = max(OD_values);


[lag_time, max_spec_growth_rate, median_od_max, delta_OD_max, doubling_time, goodness, note] = FindRegressionCurve(OD_values,time_interval, incubation_time, model, double_hump, threshold);

%report both max OD and median filtered max OD to excel 

lag_time_str = sprintf('lag time = %f, max growth = %f, doubling time = %f, Delta OD Max = %f', lag_time, max_spec_growth_rate, doubling_time, delta_OD_max);

legend(lag_time_str, 'location', 'SouthOutside');

end

function [lag_time, msgr, median_od_max, delta_OD_max, doubling_time, goodness, note] = FindRegressionCurve(OD_values, time_interval, incubation_time, model, double_hump, threshold)
%FindRegressionCurve - using the parameters specified, attempts to fit a 
% regression best fit line of the specified model to the data supplied.

msgr = 0;
note = '';
lag_time = 0;


timepoints = (0:size(OD_values)-1) * time_interval + incubation_time;


timepoints_med = timepoints;

OD_values_med = OD_values;

filter_width = 1;
for i=1+filter_width:size(OD_values)-filter_width
    sel = i - filter_width : i + filter_width;
    OD_values_med(i) = median(OD_values(sel));  
end;

min_od = min(OD_values_med);

   
[median_od_max, max_index] = max(OD_values_med);

double_hump_found = 0;
if (double_hump && max_index > 4)
    smoothing_range = 5;
    peaks = '';
    while (smoothing_range > 1 && length(peaks) < 1)
        smoothed_od_values =  smooth(OD_values, smoothing_range);
        [peaks, locs] = findpeaks(smoothed_od_values(1:min(max_index-1, length(smoothed_od_values))), 'MINPEAKDISTANCE', 3);
        j = 1;
        if (~isempty(peaks))
            peaks_thresholded = [];
            locs_thresholded = [];
            for i = 1:length(peaks)
                if (peaks(i) < median_od_max * .85)
                    peaks_thresholded(j) = peaks(i);
                    locs_thresholded(j) = locs(i);
                    j = j + 1;
                end
            end
            peaks = peaks_thresholded;
            locs = locs_thresholded;
        end        
        smoothing_range = smoothing_range - 1;
    end
    if  (length(peaks)>=1)

        [min_value, location_min] = min( smoothed_od_values(locs(length(locs)):max_index));
        location_min = location_min + locs(length(locs));

        OD_values_dh(1) = OD_values(1);
        timepoints_dh(1) = timepoints(1);

        OD_values_dh(2: length(OD_values)- location_min + 2) = OD_values(location_min: length(OD_values));
        timepoints_dh(2: length(OD_values)- location_min + 2) = timepoints(location_min: length(OD_values)); 

        OD_values = OD_values_dh;
        timepoints = timepoints_dh;
        max_index = max_index - location_min + 2;
        double_hump_found = 1;
    end
end
timepoints_orig = timepoints;
OD_values_adj(size(OD_values)) = 0;


for i=1:length(OD_values)
    if (i < max_index)
        OD_values_adj(i) = OD_values(i); % should we use the median filtered data or just the raw data? 
    else
        OD_values_adj(i) = median_od_max;      
    end
end;

max_od_adj = max(OD_values_adj);

%Appl. Environ. Microbiol. June 1990 vol. 56 no. 6 1875-1881

%and...

%0 Journal Article
%D 2000
%@ 0178-515X
%J Bioprocess Engineering
%V 23
%N 6
%R 10.1007/s004490000209
%T Development of mathematical models (Logistic, Gompertz and Richards models) describing the growth pattern of Pseudomonas putida (NICM 2174)
%U http://dx.doi.org/10.1007/s004490000209
%I Springer-Verlag
%8 2000-12-01
%A Annadurai, G.
%A Rajesh Babu, S.
%A Srinivasamoorthy, V. R.
%P 607-612
%G English

if (strcmpi(model, 'gompertz'))
      
    %Gompertz curve
    func = fittype('A * exp( -exp(-C*(x-B)))+D');

    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0 0 -1], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

    inflection_point = reg_curve.B;
    msgr = reg_curve.C;

    offset = .01;

    %is it ok to offset these values to make the log work properly?
    lag_time = inflection_point * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);


elseif (strcmpi(model, 'modgompertz'))
    %modified Gompertz curve
    func = fittype('A * exp(-exp(((B * exp(1))/ A) * (C - x) + 1)) + D');

    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0.000001 0 -1], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

    lag_time = reg_curve.B;
    msgr = reg_curve.C;


elseif (strcmpi(model, 'logistic'))
    %Logistic curve
    func = fittype('A / (1+exp(-C*(x-B))) + D');
    
    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0 0 -1], 'Upper', [100 100 2 3], 'StartPoint', [0.5 1.5 0.2 0.1]);

    inflection_point = reg_curve.B;
    msgr = reg_curve.C;

    offset = .01;

    %is it ok to offset these values to make the log work properly?
    lag_time = inflection_point * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);


    
elseif (strcmpi(model, 'modlogistic'))

    %Modified Logistic curve
    func = fittype('A / (1 + exp(((4*C)/A) * (B - x) + 2)) + D');
    
    [reg_curve, goodness] = fit(timepoints', OD_values_adj', func, 'Lower', [0 0.000001 0 -1], 'Upper', [100 100 3 3], 'StartPoint', [2 1.5 0.2 0.1]);
    size_tmp = size(OD_values_adj);
   
    lag_time = reg_curve.B;
    msgr = reg_curve.C;

%weibull does not seem to work well so disabling it
%elseif (strcmpi(model, 'weibull'))
%    %weibull
%    func = fittype('B * exp(C * (log(x) - log(A)))');
%    
%    reg_curve = fit(timepoints', OD_values_adj', func, 'Lower', [0 0 0], 'Upper', [100 100 10], 'StartPoint', [0.5 1.5 0.5]);
%    
%    inflection_point = reg_curve.B;
%    msgr = reg_curve.C;
%
%    offset = .01;
%
%    %is it ok to offset these values to make the log work properly?
%    lag_time = inflection_point * time_interval - (log((reg_curve(inflection_point) + offset) - log(OD_values_adj(1) + offset)) / msgr);
%
end;

delta_OD_max = median_od_max - min(OD_values(1:4)); %max minus initial (min of first 4 elements to remove ouliers)

%use the max specific growth rate to calc the doubling time t_d = ln(2)/u
doubling_time = log(2) / msgr;


lag_time = max([lag_time 0]);


if (goodness.rsquare < 0.98)
   note = 'Bad r^2, Check Regression Plot'; 
end

noreg = 0;


% no longer checking the absolute threshold, only delta OD (relative to
% min OD)
if (delta_OD_max < threshold)
    note = 'No Growth Detected, Check Plot';    
    lag_time = ' ';
    msgr = ' ';
    doubling_time = ' ';
    noreg = 1;    
    goodness.rsquare = ' ';
    goodness.sse = ' ';
    goodness.rmse = ' ';
    goodness.dfe = ' ';
    goodness.rsquare_adj = ' ';
end


if (double_hump==1 && double_hump_found == 0)
   note = 'No Double Hump Detected'; 
end

xlist = [lag_time lag_time];
ylist = [ -10 10];

plot(timepoints_orig, OD_values,'.')
hold on 
if (noreg == 0)
    plot(reg_curve);
end
hold on 
plot(timepoints_med, OD_values_med);
hold on 
%plot (timepoints, log_OD_values);
if (noreg == 0)
    line(xlist, ylist, 'Color', 'm', 'LineStyle', '--');
end
ylim([min(-0.2, min_od) max(1.6, 1.05*median_od_max)])
xlabel('Time (Hours)');
ylabel('OD (600nm)');

hold off
end



