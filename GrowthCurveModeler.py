import os
import glob
import xlrd
import xlsxwriter
import numpy as np
import scipy.optimize
import scipy.stats
import pylab as pl
import pdb
import math
import matplotlib.pyplot as mpl
import matplotlib.legend as mpll
import warnings


def gompertz(x, A, B, C, D):
    return A * np.exp( - np.exp( -C * (x-B))) + D



def modgompertz(x, A, B, C, D):
    return A * np.exp(-np.exp(((B * np.e)/A) * (C-x) + 1)) + D


def logistic(x, A, B, C, D):
    return A / (1 + np.exp(-C * (x - B))) + D


def modlogistic(x,A,B,C,D):
    return A / (1 + np.exp(((4 * C)/A) * (B-x) + 2)) + D


def plot_results(full_plot, full_filename, lag_time=None, timepoints=None, time = None, regression=None, OD_values=None,
                 msgr=None, doubling_time=None, delta_OD_max=None, min_od=None, median_od_max=None, possible_growth=0,
                 plot_title=''):
    fig, ax = mpl.subplots();
            
    mpl.xlabel('Time (Hours)')
    mpl.ylabel('OD (600nm)')

    ax.plot(timepoints, OD_values, 'r.')
    mpl.ylim((min(-0.2, min_od * 1.1), max(1.6, median_od_max * 1.1)))

    if (full_plot):
        lagtimestart = [lag_time, -10]
        lagtimestop = [lag_time, 10]
        #plot timepoints
        ax.plot(time, regression, 'b-')
        ax.plot(*zip(lagtimestart, lagtimestop), color='green')

        if possible_growth:
            lag_time_str = 'Possible Growth...\n'
        else:
            lag_time_str = ''
        lag_time_str = lag_time_str + "lag_time = %f\nmax growth = %f\ndoubling time = %f\nDelta OD Max = %f" % (lag_time, msgr, doubling_time, delta_OD_max)

    else:
        lag_time_str = "No Growth\nDelta OD Max = %f" % delta_OD_max

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    if median_od_max > 0.8:
        v_align = 'bottom'
        v_placement = 0.05
    else:
        v_align = 'top'
        v_placement = 0.95
    ax.text(0.95, v_placement, lag_time_str, transform=ax.transAxes, fontsize = 14, verticalalignment=v_align, horizontalalignment='right', bbox=props)
    fig.suptitle(plot_title, fontsize=16)
    mpl.savefig(full_filename + ".png", dpi =600, format="png")


    mpl.close(fig)
    #mpl.show();    
    return
    
def FindRegressionCurve(OD_values, time_interval, incubation_time, model, double_hump, threshold, full_filename, data_min, ignore_pre_min, plot_title):
    """
    FindRegressionCurve - using the parameters specified, attempts to fit a 
    regression best fit line of the specified model to the data supplied.
    """

    msgr = 0
    note = ""
    lag_time = 0
    goodness = 1
    OD_values_med = np.empty([len(OD_values)])
    timepoints = []
    for i in range(len(OD_values)):
        timepoints.append(i * time_interval + incubation_time)

#    for i in range(len(OD_values)):
#        OD_values_med[i] = max(0, OD_values[i]);
        

    timepoints_med = timepoints
     
    filter_width = 3
    
    OD_values_med = medfilt(OD_values, filter_width)

    min_od = np.min(OD_values)

    median_od_max = np.max(OD_values_med)

    max_index = np.argmax(OD_values_med)
    if max_index > 0:
        median_min_od = np.min(OD_values_med[0:max_index])
    else:
        median_min_od = OD_values_med[0]
    min_index = np.argmin(OD_values_med)



    double_hump_found = False;
    
    if (double_hump and max_index > 3):

        peaks_indices = scipy.signal.find_peaks_cwt(OD_values, np.arange(1,10))
        peaks_thresholded = []
        locs_thresholded = []
        
        if len(peaks_indices)>0:
            for i in range(len(peaks_indices)):
                if (OD_values[peaks_indices[i]] < median_od_max * 0.85):
                    peaks_thresholded.append(OD_values[peaks_indices[i]])
                    locs_thresholded.append(peaks_indices[i])

        if (len(peaks_thresholded) > 0):
            
            min_value = OD_values[peaks_thresholded[len(peaks_thresholded)-1]]
            for ele in OD_values[range(peaks_thresholded[len(peaks_thresholded)-1], len(OD_values) - 1 )]:
                min_value = min(ele, min_value)
            
            location_min = OD_values.index(min_value)
            
            OD_values_dh[0] = OD_values[0]
            timepoints_dh[0] = timepoints[0]
            
            for i in range(1, len(OD_values)-location_min +1): 
                OD_values_dh[i] = OD_values[i + location_min - 1]
                timepoints_dh[i] = timepoints[i + location_min -1]
            
            OD_values = OD_values_dh
            timepoints = timepoints_dh
            max_index = max_index - location_min + 2
            double_hump_found = True
    
    timepoints_orig = timepoints

    initial_index = 0
    if (ignore_pre_min):
        initial_index = min_index

    OD_values_adj = []
    if max_index < 4 or max_index > len(OD_values)*0.93:
        max_index = len(OD_values) + 1
        adj_vals_length = len(OD_values)
    else:
        adj_vals_length = max(max_index + 10,   math.ceil(max_index * 1.3))

    timepoints_adj = []
    for i in range(adj_vals_length):
        timepoints_adj.append(i * time_interval + incubation_time)
        if (i >= initial_index):
           if (i <= max_index):
               OD_values_adj.append(OD_values_med[i])
           else:
               OD_values_adj.append(median_od_max)
    
    
    
    
    max_od_adj = max(OD_values_adj)
    
    
    """
    Appl. Environ. Microbiol. June 1990 vol. 56 no. 6 1875-1881

    and...

    0 Journal Article
    D 2000
    @ 0178-515X
    J Bioprocess Engineering
    V 23
    N 6
    R 10.1007/s004490000209
    T Development of mathematical models (Logistic, Gompertz and Richards models) describing the growth pattern of Pseudomonas putida (NICM 2174)
    U http://dx.doi.org/10.1007/s004490000209
    I Springer-Verlag
    8 2000-12-01
    A Annadurai, G.
    A Rajesh Babu, S.
    A Srinivasamoorthy, V. R.
    P 607-612
    G English
    """
    min_d = -1
    if (data_min):
        min_d = min(0.0999,max(min(OD_values_adj), -1))

    if (model == 'gompertz'):
        fitfunc_range = gompertz
        fitfunc = gompertz
        initial_guess =  [0.5, 1.5, 0.2, 0.1]
        boundvals = ([0,0,0,min_d],[100,100,2,3])
    elif (model == 'modgompertz'):
        fitfunc_range = modgompertz
        fitfunc = modgompertz
        initial_guess = [0.5, 1.5, 0.2, 0.1]
        boundvals = ([0,0.000001,0,min_d],[100,100,2,3])
    elif (model == 'logistic'):
        fitfunc_range = logistic
        fitfunc = logistic
        initial_guess = [0.5, 1.5, 0.2, 0.1]
        boundvals = ([0,0,0,min_d],[100,100,2,3])
    elif (model == 'modlogistic'):
        fitfunc_range = modlogistic
        fitfunc = modlogistic
        initial_guess = [2, 1.5, 0.2, 0.1]
        boundvals =([0,0.000001,0,min_d],[100,100,3,3])
    else:
        print("Unsupported Model")
        return -1

    delta_OD_max = median_od_max - median_min_od
    time = pl.arange(0, (len(OD_values_adj)) * time_interval + incubation_time, 0.01)

    if (delta_OD_max < threshold * 0.7):
        note = 'No Growth Detected, Check Plot'
        lag_time = ''
        msgr = ''
        doubling_time = ''
        rsquared = ''
        rmse = ''
        plot_results(False, full_filename, 0, timepoints, time, None, OD_values, msgr, doubling_time,
                     delta_OD_max, min_od, median_od_max, plot_title=plot_title)
        return (lag_time, msgr, median_od_max, median_min_od,  delta_OD_max, doubling_time, rsquared, rmse, note)
    possible_growth_flag = 0
    if (delta_OD_max < threshold * 1.3):
        note = 'Possible Growth, Check Plot'
        possible_growth_flag = 1


    try:
        (coef, est_cov) = scipy.optimize.curve_fit(fitfunc_range, timepoints_adj[initial_index:len(OD_values_adj)], OD_values_adj, initial_guess, bounds = boundvals)
    except scipy.optimize.OptimizeWarning:
        print ('Maxed out cant estimate covariance, cannot fit curve')
        note = 'failed to fit curve'
        lag_time = ''
        msgr = ''
        doubling_time = ''
        noreg = 1
        rsquared = ''
        rmse = ''
        delta_OD_max = ''
        plot_results(False, full_filename, plot_title=plot_title)
        return (lag_time, msgr, median_od_max, median_min_od, delta_OD_max, doubling_time, rsquared, rmse, note)
    except RuntimeError:        
        print ('Maxed out calls, cannot fit curve')
        note = 'failed to fit curve'
        lag_time =''
        msgr = ''
        doubling_time = ''
        noreg = 1
        rsquared = ''
        rmse = ''
        delta_OD_max = ''
        plot_results(False, full_filename, plot_title=plot_title)
        return (lag_time, msgr, median_od_max, median_min_od, delta_OD_max, doubling_time, rsquared, rmse, note)


    #need to re-evaluate use of R^2 value.  With non-linear regression, this value is dubious.  Should perhaps use standard error in the units of OD or a pseudo R^2 value
    #(slope, intercept, rsquared, pvalue, stderr) = scipy.stats.linregress(OD_values_adj, fit_values)
    SSE = 0
    SST = 0
    adj_mean = np.mean(OD_values_adj)


    for i in range(len(OD_values_adj)):
        resid = OD_values_adj[i] - fitfunc(timepoints_adj[i], coef[0], coef[1], coef[2], coef[3])
        SST_sub = OD_values_adj[i] - adj_mean
        SSE = SSE + resid * resid
        SST = SST + SST_sub * SST_sub

    rsquared = 1- SSE/SST
    rmse = np.sqrt( SSE / len(OD_values_adj))


    if (model == 'gompertz' or model == 'logistic'):
        inflection_point = coef[1]
        msgr = coef[2]
        #offset = 1
        #lag_time = inflection_point * time_interval - (np.log(fitfunc(inflection_point,coef[0],coef[1],coef[2], coef[3]) + offset) - np.log(OD_values_adj[0] + offset)) / msgr
        lag_time = inflection_point * time_interval - (fitfunc(inflection_point,coef[0],coef[1],coef[2], coef[3]) - fitfunc(0,coef[0],coef[1],coef[2],coef[3])) / msgr

    elif (model == 'modgompertz' or model == 'modlogistic'):
        lag_time = coef[1]
        msgr = coef[2]


        
    doubling_time = np.log(2) / msgr
    
    lag_time = max(lag_time, 0)
    
    
 
    
    if (double_hump==1 and double_hump_found == 0):
        note = 'No Double Hump Detected'
        

    regression = fitfunc(time, coef[0],coef[1], coef[2], coef[3])

    plot_results(True, full_filename, lag_time, timepoints, time, regression, OD_values, msgr, doubling_time, delta_OD_max, min_od, median_od_max, possible_growth_flag, plot_title=plot_title)
    
    #plot curve
    
    #plot line on lag_time
    
    #scale, and label
    
    return (lag_time, msgr, median_od_max, median_min_od, delta_OD_max, doubling_time, rsquared, rmse, note);
    

def medfilt (x, k):
    """Apply a length-k median filter to a 1D array x.
    Boundaries are extended by repeating endpoints.
    """
    assert k % 2 == 1, "Median filter length must be odd."
    assert x.ndim == 1, "Input must be one-dimensional."
    k2 = (k - 1) // 2
    y = np.zeros ((len (x), k), dtype=x.dtype)
    y[:,k2] = x
    for i in range (k2):
        j = k2 - i
        y[j:,i] = x[:-j]
        y[:j,i] = x[0]
        y[:-j,-(i+1)] = x[j:]
        y[-j:,-(i+1)] = x[-1]
    return np.median (y, axis=1)
    
    

def MicrobialKinetics(OD_values, time_interval, incubation_time, threshold, model, double_hump, full_filename, data_min, ignore_pre_min, plot_title):
    """
    MicrobialKinetics -  For a specific dataset, attempts to determine a group
    of statistics on the growth which occurred.   Plots the dataset along the 
    way, with the best fit regression.
    """
    # time interval assumed to be half hour time blocks


    #set to zero initially   
    max_od  = np.max(OD_values)
    max_location = np.argmax(OD_values)
    min_od = np.min(OD_values)

    (lag_time, max_spec_growth_rate,median_od_max, median_od_min, delta_OD_max, doubling_time, rsquared, rmse, note) = FindRegressionCurve(
                    OD_values, time_interval, incubation_time, model, double_hump, threshold, full_filename, data_min, ignore_pre_min, plot_title)
    
    
    #report both max OD And median filtered max OF to excel


    
    
    #legend(lag_time_str, 'location', 'SouthOutside');
    
    
    return (lag_time, max_spec_growth_rate, max_od, min_od, median_od_max, median_od_min, delta_OD_max, doubling_time, rsquared,rmse, note)
    
    
    
def GrowthCurveModeler( file_or_dir, **varargin):
    """
    GrowthCurveModeler - Calculates some metrics for growth curves, as well as graphing
        a regression curve of the data points. 

    Required parameter: 
    
    file_or_dir - input data file or directory of input files, must be formatted properly

    Optional parameters:

    MaxTimepoint - value {default is total specified in dataset}
      Final timepoint of the dataset 

    Threshold - value {default 0.3}
      Threshold which describes the minimum OD reading to signify growth (default 0.3)
 
     Model - [{'modlogistic'} | 'gompertz' | 'logistic' | 'modgompertz']
       Regression model to plot the curves

     PreIncubationTime - value {default 1.0}
       If you incubate your cells before recording timepoints,
       this will appropriately shift your data points in time
       for lag time calculations 

     DoubleHump - [True | {False}] 
       Flag to indicate double hump processing, this expects all datasets to include a
       double/multi hump and should remove all growth humps before the "main curve" 
 
     DataMin - [ True | {False}]
        Data regression minimum is min data point

     IgnorePreMin - [ True | {False}]
        Data regression ignores data before the minimum data point (prior to the maximum datapoint)

     RSquaredFlag - value {default .97}
        Cutoff for "Good" regression fit

      Examples:
          GrowthCurveModeler('dataset.xlsx', DoubleHump=True, Threshold=0.2);

         GrowthCurveModeler('.', IncubationTime=1.5);
 
         GrowthCurveModeler('folder_containing_xlsxfiles');
    """
    print(file_or_dir)
    if (os.path.isdir(file_or_dir)):
        runnable_files = glob.glob(file_or_dir + "/*.xlsx")
        for r in runnable_files:
            GrowthCurveModeler(r, **varargin)
        return

    max_timepoint = -1
    growth_threshold = 0.3
    model = 'modlogistic'
    incubation_time = 1.0
    double_hump = 0
    data_min = 0
    ignore_pre_min = 0
    r2_good_fit_cutoff = 0.97

    for k,v in varargin.items():
        if (k=='MaxTimepoint'):
            max_timepoint = v
        if (k=='Threshold'):
            threshold = v
        if (k=='Model'):
            model = v
        if (k=='PreIncubationTime'):
            incubation_time = v
        if (k=='DoubleHump'):
            double_hump = v
        if (k=='DataMin'):
            data_min = v
        if (k=='IgnorePreMin'):
            ignore_pre_min = v
        if (k=='RSquaredFlag'):
            r2_good_fit_cutoff = v
    
    (path, file) = os.path.split(file_or_dir)
    (stub, ext) = os.path.splitext(file)
    
    if (path):
        path = path + "/"
        
    plots_folder = path + "results/"
        
    if (not os.path.exists(plots_folder)):
        os.mkdir(plots_folder)
        
    plots_folder = plots_folder + stub + " plots/"

    if (not os.path.exists(plots_folder)):
        os.mkdir(plots_folder)

    
    workbook = xlrd.open_workbook(file_or_dir)
    sheet = workbook.sheet_by_index(0)
    num_columns = sheet.ncols
    
    output = ('Sugar', 'Strain', 'Lag Time (hours)', 'Max Specific Growth Rate (1/hours)',\
            'Doubling Time (hours)', 'Max OD', 'Max OD (Median Filtered Data)', 'Min OD', 'Min OD (Median Filtered Data)', 'Delta OD (Median Filtered Data)', 'Notes', 'R^2', 'RMSE')
    
    time_interval = 0.5
    sugar = ''
    strain_count = 0
    sugar_count = 0

    first_sugar = True
    Sugars = []
    Start_idxs = []
    Strain_counts = []
    lag_times = []

    output_file = path + 'results/' + stub + ' results.xlsx'
    output_workbook = xlsxwriter.Workbook(output_file)
    output_sheet = output_workbook.add_worksheet("Results")
    green_fill = output_workbook.add_format()
    green_fill.set_bg_color('#80D040')
    red_fill = output_workbook.add_format()
    red_fill.set_bg_color('#FF5050')
    red_fill.set_align('center')
    no_fill = output_workbook.add_format()
    no_fill.set_align('center')
    red_text = output_workbook.add_format()
    red_text.set_color('red')
    red_text.set_align('center')

    output_sheet.conditional_format('L2:L600', {'type': 'cell', 'criteria': 'between',  'maximum': r2_good_fit_cutoff,'minimum': 0.0000001,  'format': red_fill})
    output_sheet.conditional_format('L2:L600', {'type': '2_color_scale','min_color': "#FFA550",'max_color': "#80D040", 'min_type': 'num', 'max_type':'num','min_value': r2_good_fit_cutoff,'max_value':1.0})
#   output_sheet.conditional_format('L2:L600', {'type': 'cell','criteria': '>', 'value': r2_good_fit_cutoff,'minimum': 0.0000001,  'format': green_fill})

    output_sheet.conditional_format('K2:K600', {'type': 'no_blanks', 'format': red_text})
    output_sheet.freeze_panes(1,2)
    header_format = output_workbook.add_format()
    header_format.set_text_wrap(1)
    header_format.set_bold(1)
    header_format.set_align('center')
    data_format = output_workbook.add_format()
    data_format.set_align('center')

    output_sheet.write_row(0,0, output, header_format)
        
    title_data0 = sheet.row(0)
    title_data1 = sheet.row(1)
    for i in range(num_columns):
        data_column = sheet.col_slice(i, 2)
        data_column_values = []
        for d in data_column:
            data_column_values.append(float(d.value))
        if (sheet.cell(0,i)):
            if (not first_sugar):
                Strain_counts.append(strain_count)
            strain_count = 0
            sugar_start_idx = i
            first_sugar = False
            sugar_count = sugar_count + 1
            sugar = title_data0[i].value
            Start_idxs.append(sugar_start_idx)
            Sugars.append(sugar)
        if (title_data1[i].value=="Time"):
            time_interval = data_column_values[1] - data_column_values[0]
        else:
            strain_count = strain_count + 1
            strain = title_data1[i].value
            sugar_folder = plots_folder + "/" + sugar
            if (not os.path.exists(sugar_folder)):
                os.mkdir(sugar_folder)
            full_filename = sugar_folder + "/" +sugar + '-' + strain
            
            if (max_timepoint < 0):
                (lagtime, max_u, OD_max, OD_min, median_OD_max, median_OD_min, delta_OD_max, doubling_time, rsquared, rmse, note) = MicrobialKinetics(
                    np.array(data_column_values, dtype=float), time_interval, incubation_time, growth_threshold, model,
                    double_hump, full_filename, data_min, ignore_pre_min, sugar + '-' + strain)
            else:
                (lagtime, max_u, OD_max, OD_min, median_OD_max, median_OD_min, delta_OD_max, doubling_time, rsquared, rmse, note) = MicrobialKinetics(
                    np.array(data_column_values[0:max_timepoint / time_interval]), time_interval, incubation_time,
                    growth_threshold, model, double_hump, full_filename, data_min, ignore_pre_min, sugar + '-' + strain)
            lag_times.append(lagtime)

            if (not isinstance(rsquared, str) and rsquared  < r2_good_fit_cutoff):
                if (len(note) > 0):
                    note = note + '; '
                note = note + 'Poor Regression'

            output_sheet.write_row(i, 0, (sugar,strain,lagtime, max_u, doubling_time, OD_max, median_OD_max, OD_min, median_OD_min, delta_OD_max, \
                                          note, rsquared, rmse), data_format)

            
            
            #plot title = name         
             
            #save plot to sugar folder
            #close plot
    
    
    
 

    


    