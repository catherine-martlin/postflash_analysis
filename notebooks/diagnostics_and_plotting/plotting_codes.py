from astropy.time import Time
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, Range1d
from bokeh.layouts import row, column, gridplot
from bokeh.models.widgets import Tabs, Panel
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

def plot_mean_normalized_and_linear_fit(pandas_df, legend_label):
    # normalized to the mean of first observation
    x_datetime = pd.to_datetime(pandas_df['obsdate'] + ' ' + pandas_df['obstime'])
    sorted_x = x_datetime.sort_values()
    first_rootname = sorted_x.index[0]
    x = Time(x_datetime,format='datetime64')
    x.format = 'mjd'
    x = x.to_value('mjd', 'float')

    mean_mean = np.mean(pandas_df['norm_by_dur'])
    y = pandas_df['norm_by_dur']/pandas_df['norm_by_dur'][first_rootname]

    # Make and fit a linear regression model
    def line(x, a, b):
        return a * x + b

    #x = x.reshape(-1, 1)
    popt, pcov = curve_fit(line, x, list(y))
    print("Slope: ", popt[0])
    perr = np.sqrt(np.diag(pcov))
    print(perr)
    print(perr[0]/popt[0])

    # Find the slope and intercept from the model
    slope = popt[0] # Takes the first element of the array
    intercept = popt[1]

    # Make the regression line
    y_pred = []
    for i in x:
        i = float(i)
        temp_y = slope*i + intercept
        y_pred.append(temp_y)

    p = figure(title = 'Mean PF value Normalized to Mean First Observation and Duration',x_axis_type='datetime',plot_width=800, plot_height=300)
    if 'Full-frame' in legend_label:
        p.y_range=Range1d(0.95, 1.05)
    else:
        p.y_range=Range1d(0.95, 1.05)
    p.circle(x_datetime, y,size = 4, legend_label = legend_label)
    p.line(x_datetime,y_pred,color='red')#,legend_label='y = ' + str(slope) +' x + '+ str(round(intercept,6)))
    p.xaxis.axis_label = "Observation Date"
    p.yaxis.axis_label = "Normalized Mean Post-flash Value"
    show(p)
