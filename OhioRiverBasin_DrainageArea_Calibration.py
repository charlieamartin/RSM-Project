'''
This program reads an excel spreadsheet containing Ohio Basin USGS gages.....
'''
import csv
import urllib.request, urllib.error, urllib.parse
import re
import matplotlib.pyplot as plt
import numpy as np

from math import log10, pow
from OhioRiverBasin_DrainageArea_Functions import compPowerLaw, biasCorrect, isFloat, PlotHistogram



# Main Code
#
# you must provide the following information
#
filename = 'newCalibrationSites.csv'

#
# this is the main code
#

sitenum = []
drainarea = []
plaws = []

# Drainage Areas Smaller than 200
plawl200 = []
sedxl200 = []
sedyl200 = []

# Drainage Areas Greater than 200
plawg200 = []
sedxg200 = []
sedyg200 = []
#error_all = []
#rela_error_all = []

#
# open and read the site locations from a csv file
#
with open(filename) as csvfile:

    readCSV = csv.reader(csvfile, delimiter=',')

    for row in readCSV:
       sitenum.append(row[0])
       drainarea.append(row[6])


file_csv1 = open('Global_Error_Calibration.csv', 'w')
file_csv1.write("{0},{1}\n".format("Absolute Error", "Relative Error"))
file_csv2 = open('lessthan200_Error_Calibration.csv', 'w')
file_csv2.write("{0},{1}\n".format("Absolute Error", "Relative Error"))
file_csv3 = open('greaterthan200_Error_Calibration.csv', 'w')
file_csv3.write("{0},{1}\n".format("Absolute Error", "Relative Error"))

#
#  get data from USGS corresponding to sites in the csv file
#
#for i in range(1,len(sitenum)):
    #siteno = sitenum[i]
    #drainsize = drainarea[i]

for siteno, drainsize in zip(sitenum[1:], drainarea[1:]):
    if not isFloat(drainsize):
        print('   ')
        print(siteno, 'with drainange area = ', drainsize)
        print('   ')
        continue

    drainsize = float(drainsize)
    sedx = []
    sedy = []
    Qhat = []
    sedqx = []

    #Field measurements
    '''
    url='https://waterdata.usgs.gov/nwis/measurements?site_no='
    url=url + siteno
    url=url + '&agency_cd=USGS&format=rdb_expanded'
    page = urllib2.urlopen(url)
    page_content = page.read()
    ifile=open('field_'+siteno+'.txt', 'w')
    ifile.write(page_content)
    ifile.close
    '''

    #Sediment measurements
    codes = 'parm_cds=00061%2C70331%2C80154&qw'   #retrieving 00061, 70331, and 80154

    url = 'https://nwis.waterdata.usgs.gov/nwis/qwdata/?site_no='
    url = url + siteno
    url = url + '&agency_cd=USGS&inventory_output=0&rdb_inventory_output=file&TZoutput=0&pm_cd_compare=Greater%20than&radio_parm_cds=parm_cd_list&radio_multiple_'
    url = url + codes
    url = url + '_attributes=0&format=rdb&qw_sample_wide=wide&rdb_qw_attributes=0&date_format=YYYY-MM-DD&rdb_compression=value&submitted_form=brief_list'
    page = urllib.request.urlopen(url)

    # Decode required to ensure the page can be read
    page_content = page.read().decode('utf-8')

    #
    # look for discharge and suspended sediment flags
    #
    p6flag = False
    p7flag = False
    p8flag = False
    i6 = 0
    i7 = 0
    i8 = 0
    llist = page_content.split('\n')
    for line in llist:
        # split the line into words separated by tabs
        words = line.replace('\r', '').split('\t')

        # discharge
        if re.match("#\s+P00061", line):
            p6line = line
            p6flag = True
        # suspended sediment
        elif re.match("#\s+P70331", line):
            p7line = line
            p7flag = True
        # suspended sediment concentration
        elif re.match("#\s+P80154", line):
            p8line = line
            p8flag = True

        if p6flag and p7flag and p8flag:
            print(siteno)
            print(p6line)
            print(p7line)
            print(p8line)
            #print("{0:>16}{1:>13}{2:>11}".format('P00061','P70331','P80154'))
            break

    #
    # if the data is not there, go to the next file
    #
    if not (p6flag and p7flag and p8flag):
        continue

    #
    # go through the lines again to pull out the data
    #
    for line in llist:
        # split the line into words separated by tabs
        words = line.replace('\r', '').split('\t')

        #
        # find the header record
        #
        if p6flag and p7flag and p8flag and re.match("agency_cd", line):

            # find the indices of the discharge and sediment data
            i6 = words.index('p00061')
            i7 = words.index('p70331')
            i8 = words.index('p80154')

        # the data lines start with USGS
        if re.match("\AUSGS", line):
            if len(words[i6]) > 0 and len(words[i8]) > 0:
                w6 = words[i6].strip('EA<>')
                w7 = words[i7].strip('EA<>')
                w8 = words[i8].strip('EA<>')
                dis = float(w6)
                con = float(w8)
                #print(" sed {0:12.3f}{1:>12}{2:12.3f}".format(dis,"  -  ",con))
                if (dis != 0.0) and (con != 0.0):
                   sedx.append(dis)
                   sedy.append(con)
    #
    # perform linear regression by polynomial fitting
    #

    a, b, pts = compPowerLaw(sedx, sedy, 5)
    plaws.append(pts)
    if drainsize >= 200:
        print('> 200 drainsize = ', drainsize)
        plawg200.append(pts)
        sedxg200.append(sedx)
        sedyg200.append(sedy)
    else :
        print('< 200 drainsize = ', drainsize)
        plawl200.append(pts)
        sedxl200.append(sedx)
        sedyl200.append(sedy)

    print('      1st power law coeffs a, b = ', a, b)


plawxall = []
plawyall = []
for plaw in plaws:
    plawxall.extend(plaw[0])
    plawyall.extend(plaw[1])

a, b, pts = compPowerLaw(plawxall, plawyall, 5)
print('   ')
print('Universal power law coeffs a, b = ', a, b)

#plt.figure()
#plt.loglog(plawxall, plawyall, 'g*')
#plt.loglog(pts[0], pts[1], 'g*-')

qhat = biasCorrect(plawxall, plawyall, a, b)

#Compute Universal Errors
np_plawxall = np.array(plawxall)
np_plawyall = np.array(plawyall)

#convert from (mg/L)*(ft^3/sec) to (ton/day) multiply by 0.002698
error_all = 0.002698 * (qhat * np_plawxall - np_plawyall * np_plawxall)
rela_error_all = error_all / (np_plawyall * np_plawxall * 0.002698)

for i in range(len(error_all)):
    file_csv1.write("{0},{1}\n".format(error_all[i],rela_error_all[i]))

PlotHistogram(error_all, 'Histogram of Error for Global Data')

#a,b, pts = compPowerLaw(plawxall,qhat,5)
#plt.loglog(pts[0], pts[1], 'k*-')

# apply the bias correction
plawlxall = []
plawlyall = []
for plawl in plawl200:
    plawlxall.extend(plawl[0])
    plawlyall.extend(plawl[1])

a, b, pts = compPowerLaw(plawlxall, plawlyall, 5)
print('   ')
print('    < 200 power law coeffs a, b = ', a, b)

#plt.loglog(plawlxall, plawlyall, 'b*')
#plt.loglog(pts[0], pts[1], 'b*-')

qhat = biasCorrect(plawlxall, plawlyall, a, b)

#a,b, pts = compPowerLaw(plawlxall,qhat,5)
#plt.loglog(pts[0], pts[1], 'c*-')

#Compute Universal Errors
np_plawlxall = np.array(plawlxall)
np_plawlyall = np.array(plawlyall)

#convert from (mg/L)*(ft^3/sec) to (ton/day) multiply by 0.002698
error_all = 0.002698 * (qhat * np_plawlxall - np_plawlyall * np_plawlxall)
rela_error_all = error_all / (np_plawlyall * np_plawlxall * 0.002698)

for i in range(len(error_all)):
    file_csv2.write("{0},{1}\n".format(error_all[i],rela_error_all[i]))

PlotHistogram(error_all, "Histogram of Error for Sites < 200 sqml")


plawgxall = []
plawgyall = []
for plawg in plawg200:
    plawgxall.extend(plawg[0])
    plawgyall.extend(plawg[1])

a, b, pts = compPowerLaw(plawgxall, plawgyall, 5)
print('   ')
print('    > 200 power law coeffs a, b = ', a, b)

#plt.loglog(plawgxall, plawgyall, 'r*')
#plt.loglog(pts[0], pts[1], 'r*-')

qhat = biasCorrect(plawgxall, plawgyall, a, b)

#a,b, pts = compPowerLaw(plawgxall,qhat,5)
#plt.loglog(pts[0], pts[1], 'm*-')

#Compute Universal Errors
np_plawgxall = np.array(plawgxall)
np_plawgyall = np.array(plawgyall)

#convert from (mg/L)*(ft^3/sec) to (ton/day) multiply by 0.002698
error_all = 0.002698 * (qhat * np_plawgxall - np_plawgyall * np_plawgxall)
rela_error_all = error_all / (np_plawgyall * np_plawgxall * 0.002698)

for i in range(len(error_all)):
    file_csv3.write("{0},{1}\n".format(error_all[i],rela_error_all[i]))

PlotHistogram(error_all, 'Histogram of Error for Sites > 200 sqml')

# Shows all plots to the user
plt.show()

file_csv1.close()
file_csv2.close()
file_csv3.close()
