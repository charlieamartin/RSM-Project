import numpy as np
import matplotlib.pyplot as plt

from math import log10


#
# compute linear regression of log10 values
#
def compPowerLaw(xin, yin, npts):
    '''
    Compute power law coefficients by computing the linear regression of log10 of input values
    Inputs:
      list: xin - x values
      list: yin - y values
      int: npts - number of points for plotting
    Output:
      tuple: (a, b, pts)
      float: a, b - coefficients for ax + b,
      list: pts - points for plotting as a list of lists [[x1,x2,...],[y1,y2,...]]
    '''
    xl = [log10(x) for x in xin]
    yl = [log10(y) for y in yin]
    a, b = np.polyfit(xl, yl, 1)

    # compute points on line for plotting - evenly spaced in log space
    xmin = min(xin)
    xmax = max(xin)
    dx = (log10(xmax) - log10(xmin)) / (npts - 1)
    plawx = [10.**(log10(xmin) + i * dx) for i in range(npts)]
    # Non Log space calculation
    #dx = (xmax- xmin)/(npts-1)
    #plawx = [xmin + i*dx for i in range(npts)]
    plawy = [(10.**b) * (x**a) for x in plawx]
    pts = [plawx, plawy]

    return (a, b, pts)


def biasCorrect(sedx, sedy, a, b):
    '''
    Calculate bias corrector as in Kao et al for a set of x,y data and the corresponding power
    law coefficients for that data.
    :rtype : list
    Inputs:
        list: sedx - x value of data points
        list: sedy - y value of data points
        float: a - power law coefficient of form ax + b
        float: b - power law coefficient of form ax + b
    Outputs:
        list: qhat - new y values based on bias correction
    '''
    np_sedx = np.array(sedx)
    np_sedy = np.array(sedy)
    tmp = (10**b) * np_sedx**a

    #compute epsilon from Eq (4) of Kao et al
    epsilon = np_sedy - tmp

    #Calculate beta bias-correction factor Eq. (5) of Kao et al
    beta = sum(epsilon) / sum(tmp)
    print("beta =", beta)

    # compute modified concentration values from the paper equation using the x values of the data points
    qhat = (1+beta) * (10.**b) * np_sedx**a

    #file_csv.write("{0},{1},{2},{3},{4},{5}\n".format(siteno, a, b, sedmin[-1],sedmax[-1],beta))

    return qhat

#
#  check if a string is a number
#
def isFloat(str):
    try:
        float(str)
        return True
    except ValueError:
        return False



def PlotHistogram(error_all, plot_title):
    plt.figure()
    n, bins, patches = plt.hist(error_all, 100, facecolor = 'green', alpha = 1.0)
    bc = 0.5 * (bins[:-1] + bins[1:])

    # Add a best fit line
    line = plt.plot(bc, n, 'r--', linewidth = 1)

    plt.xlabel('Error Range')
    plt.ylabel('Number of Measurements')
    plt.title(plot_title)
    plt.grid(True)
