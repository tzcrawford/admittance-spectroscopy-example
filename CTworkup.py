#!/usr/bin/python
#execute this with something like `./CT_workup.py > CT_workup.log 2>&1` to save results and errors to file

#Takes output data from ct.dat and works up the data using admittance spectroscopy technique to give the defect energy and the attempt to escape frequency
#need this equation: ln(f/T0^2) = -Ed/k * 1/T + ln(nu0)
#where Ed is the defect energy and nu0 is the attempt to escape frequency 

import numpy
import math
from numpy import exp, log, sqrt
    #note numpy log can only do base e
import copy
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, Data, RealData

### File Locations (EDIT THIS)

locInput="/path/to/example_of_code/admittance_spectroscopy/data/ct.dat"

#directory where output files will be saved
locOutput="/path/to/example_of_code/admittance_spectroscopy/data/CT_out.tsv"

###Definitions (EDIT THIS)
boltzmann=8.617333262145e-5 #eV/K
relative_err_in_C=0.02 #2%

#array of temperature touples where you can see a C step. For example, if you see a C step over 120 to 160K and another over 180 to 220K, you would input [(120,160),(180,220)]
stepRanges=[(100,250)] 

delimiter=str("\t")
ColumnSetTemp=1 #temperature the sample stage was set to, assuming in K
ColumnSampleTemp=2 #temperature actually measured at the sample, assuming in K
ColumnFrequency=3 #assuming in Hz
ColumnCapacitance=4
ConverstionFactorCapacitance=1 #Factor to convert capacitances in data file to nF/cm^2
ColumnConductance=5 #We don't use
NumberHeaderLines=1
#we are assuming data file has temperatures and frequencies in decending order, with every frequency at a given temp done before the next temp begins


#summary of what we are doing over the course of the whole script, when you execute the script as a whole, it will run this function using the global parameters at the top of the script
def main(locInput,locOutput,stepRanges):
    
    #data is an array of arrays: [ set temperatures , sample temperatures , frequencies , capacitances ]. Each row is a data point taken during the experiment (each line in input file).
    print("Collecting data from input file...")
    data=ReadData(locInput)
   
    Eds=list();nu0s=list(); #will hold the results for each step range

    #Remember stepRange is just a touple of which temperatures to do the calculation over (where the step occurs)
    #rangeNumber is an index number for which stepRange we are on
    for rangeNumber,stepRange in enumerate(stepRanges):

        print("\n\n-------------------------------------------------------------------")
        print("-------------------------------------------------------------------")
        print("Performing admittance spectroscopy analysis on")
        print("   temperature step range " + str(stepRanges) + " K...")
        print("-------------------------------------------------------------------")
        print("-------------------------------------------------------------------\n\n")

        #T0s contains each temperature for the inflection point at a given freqeuncy; each of these will be outputted to file
        #we also collect the other parameters for each frequency
        frequencies,T0s,err_T0s,Cmins,err_Cmins,Ls,err_Ls,ks,err_ks,chisq_T0s=getT0s(data,stepRange,rangeNumber)

        #arrheniusX & Y are the vector arrays, with each element corresponding to one data point on the arrhenius plot, ln(f/T0^2) vs 1/T0
        arrheniusX,arrheniusY,err_arrheniusX,err_arrheniusY=getarrhenius(frequencies,T0s,err_T0s)
        Ed,nu0,err_Ed,err_nu0,chisq_arrhenius=arrheniusfitting(arrheniusX,arrheniusY,err_arrheniusX,err_arrheniusY)
        

        writeData(rangeNumber,Ed,nu0,err_Ed,err_nu0,chisq_arrhenius,frequencies,T0s,err_T0s,chisq_T0s,arrheniusX,arrheniusY,err_arrheniusX,err_arrheniusY)
        
        Eds.append(Ed);nu0s.append(nu0) # We return these for each step range, although they have already been printed in respective files for each step range.
    return Eds, nu0s # We don't have to save these when we call main(), but it might be useful if calling the function otherwise.


#the function we are using to fit the capacitance vs temperature step for each frequency
def logisticfn(T,T0,Cmin,L,k):
    return L/(1 + exp(-k*(T-T0)) ) + Cmin
numberofparameters_logistic=4

#the function we are using to fit arrhenius plot
#we can't use least squares regression with error in both x and y variables, so we implement orthogonal distance regression
# https://stackoverflow.com/questions/26058792/correct-fitting-with-scipy-curve-fit-including-errors-in-x
# https://docs.scipy.org/doc/scipy/reference/odr.html
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.Model.html#scipy.odr.Model
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.RealData.html
#def linearfn(x,m,b): #want to use this format for scipy.curve_fit (least squares regression)
#    return m*(x-x0)+(b-b0)
def linearfn(beta,x): #want to use this format for scipy.odr
    #x0 and b0 assist with the fit in case the data is very far from 0, thus the fitting occurs around (x0,b0) instead of (0,0)
    return (beta[0]-b0)+beta[1]*(x-x0)
numberofparameters_linear=2



#get data from output file
def ReadData(locInput):
    
    #now reading the file
    f=open(locInput,"r")
    lines = list()
    lines = f.readlines() #saving the entire contents of the datafile as a list where each line is now an element
    f.close()

    ##now separating that data into the workable arrays
    for i in range(NumberHeaderLines):
        del lines[0] #delete the header which holds something like "Capacitance in nF/cm^2"
   
    #data is a list of lists where 
        #element 0 is set temps
        #element 1 is sample temps
        #element 2 is frequencies
        #element 3 is capacitances
    data=[ list() for i in range(4) ] #generating 4 element list

    for i in lines:
        line=i.strip('\n').split(delimiter) #separating that line by delimiter characters and deleting any newline characters
        settemp=float(line[ColumnSetTemp-1]) #should be in K
        sampletemp=float(line[ColumnSampleTemp-1]) #should be in K
        frequency=float(line[ColumnFrequency-1]) #should be in Hz
        capacitance=float(line[ColumnCapacitance-1]) * ConverstionFactorCapacitance #This should now be in nF/cm^2
        data[0].append(settemp)
        data[1].append(sampletemp)
        data[2].append(frequency)
        data[3].append(capacitance)

    return data


#function gets the temperature to use on arrhenius plot for each frequency
def getT0s(data,stepRange,rangeNumber):
   
    #remember data is a list of lists where 
        #element 0 is set temps
        #element 1 is sample temps
        #element 2 is frequencies
        #element 3 is capacitances

    #print("Preparing data for analysis on step range " + str(stepRange) + "...")

    #We need to remove data points for temperatures out of step range.
    #This will still keep the same order as data when it was inputted, but now rows are removed.
    #Going to iterate over data in the reverse direction so we can change the array while we are looping over it.
    #Note that we are changing data in this function, but does not change the version of data in main() because data is local to this function and isn't global. 
    #It is still decending in temperature and freq.
    dummy=list(range(len(data[0]))) #the list of elements numbers in data 0: [0,1,2,...]
    dummy.reverse() #example: [32,31,30,...,2,1,0]
    for i in dummy:
        #remove the row in data if the temperature is too low for stepRange
        if data[1][i] < stepRange[0]: 
            for j in data: #each column
                del j[i] #think del data[*][i] if you can regex
        
        #remove the row in data if the temperature is too high for stepRange
        if data[1][i] > stepRange[1]: 
            for j in data: #each column
                del j[i] #think del data[*][i] if you can regex
    
    frequencies=list()
    for i in data[2]:
        frequencies.append(round(i,3)) #don't want to have some kind of issue where frequency is rounded differently at some other point
    #delete duplicates
    dummy=list(range(len(frequencies)))
    dummy.reverse()
    for i in dummy:
        if frequencies[i] in (frequencies[:i] + frequencies[i+1:]): #if this frequency is somewhere else in frequencies
            del frequencies[i]
    #now frequencies is a list containing each frequency, rounded to 3 decimal digits, in descending order without any duplicates
    
    Cmin_guess=data[3][-1] #estimate for the lower bound of the sigmoid on this step
    Cmax_guess=data[3][0] #estimate for the upper bound of the sigmoid on this step


    #these will hold the parameters determined at each frequency
    T0s=list();err_T0s=list();Cmins=list();err_Cmins=list();Ls=list();err_Ls=list();ks=list();err_ks=list();chisq_T0s=list()
   
    for frequency in frequencies:

        print("\n\n-------------------------------------------------------------------")
        print("Fitting capacitance vs temperature data for frequency = " + str(frequency) + "Hz")
        print("-------------------------------------------------------------------\n")

        data_at_freq=copy.deepcopy(data) #this will be portion of data that is at this frequency
        #we cannot do just data_at_freq=data because changes in one will change the other
        dummy=list(range(len(data[0])))
        dummy.reverse()
        for i in dummy:
            #remove the row in data_at_freq if the frequency doesn't match
            if( round(data_at_freq[2][i],3) != frequency):
                for j in data_at_freq: #each column
                    del j[i] #think del data_at_freq[*][i] if you know regex
        del dummy
        x=data_at_freq[1] #sample temps at this freq
        y=data_at_freq[3] #capacitances at this freq
        
        #take the lowest capacitance, that which was taken at the lowest temperature
        C0_guess=(Cmin_guess+Cmax_guess)/2 #estimate for the capacitance that happens at T0

        T0_guess=x[min( range(len(x)), key = lambda i: abs(y[i]-C0_guess) )]
        #this command is bulky.
        #You probably aren't familiar with lambda, but it's essentially shorthand for a function that you don't want to permanently define
        #within the min() function call, we are defining a function key(i) as an optional parameter for min()
        #key works like in sorted(), where it's just a function to describe how we want to sort
        #in this case, we want to minimize the expression abs(y[i]-C0_guess) for each i in range(len(x))
        #and key would be equivalent to a function of parameter i with a return abs(y[i]-C0_guess)
        #so in the end, T0_guess is the x[i] associated with the y[i] closest to C0_guess.
        
        
        #I'm not sure the best way to estimate k, but here is my thought process
        #k is the logistic growth rate, which I think can be thought of like this:
        #The derivative of the logistic curve is something gaussian-like where the slope of C vs T is maximum at T0.
        #I think k can be thought of as the temperature distance from T0 where you are one standard deviation away
        #68% of the change in capacitance is over period of temperature one standard deviation away from T0
        #In other words, the first 16% of the change in capacitance occurs up to the temperature one standard deviation below T0
        C1sig_low=Cmin_guess + 0.16*(Cmax_guess - Cmin_guess) #capacitance that is produced 1 standard deviation below T0
        C1sig_high=Cmin_guess + (1-0.16)*(Cmax_guess - Cmin_guess)
        T1sig_low=x[min( range(len(x)), key = lambda i: abs(y[i]-C1sig_low) )] #gets the temperature closest to the one that produces C1sig_low
        T1sig_high=x[min( range(len(x)), key = lambda i: abs(y[i]-C1sig_high) )]
        k_guess=abs((Cmax_guess-Cmin_guess))/( T1sig_high - T1sig_low )

        #initialguess=[T0,Cmin,L,k]
        initialguess=[T0_guess,Cmin_guess,Cmax_guess-Cmin_guess,k_guess]
        print("Here is the initial guess array for [T0,Cmin,L,k]")
        print(str(initialguess)+"\n")
        s=[relative_err_in_C*i for i in y ] # estimated error for each capacitance in the fit

        #fitting capacitance data vs sample temperature with the logisticfn using the initial parameter guess and estimated err above
        fit=curve_fit(f=logisticfn,xdata=x,ydata=y,p0=initialguess,sigma=s,absolute_sigma=True)
        popt,pcov=fit #the fit will return two matrices. popt = the fit parameters, pcov = the covariance matrix (err in parameters)
        fit_T0,fit_Cmin,fit_L,fit_k=popt
        err_fit_T0=sqrt(numpy.diag(pcov)[0])
        err_fit_Cmin=sqrt(numpy.diag(pcov)[1])
        err_fit_L=sqrt(numpy.diag(pcov)[2])
        err_fit_k=sqrt(numpy.diag(pcov)[3])

        print("Here are final fit parameters:")
        print("T0 = " + str(fit_T0) + " +/- " + str(err_fit_T0) + " K")
        print("Cmin = " + str(fit_Cmin) + " +/- " + str(err_fit_Cmin) + " nF/cm^2")
        print("L = Cmax-Cmin = " + str(fit_L) + " +/- " + str(err_fit_L) + " nF/cm^2")
        print("k = " + str(fit_k) + " +/- " + str(err_fit_k) + " K")
        print("") #blank line
        print("Put this equation in gnuplot etc. Using x for temperature in K")
        print(str(frequency) + ":  C(x) = ( " + str(fit_L) + " ) / ( 1 + exp(-"+ str(fit_k) +"*(x - " + str(fit_T0) + "))) + " + str(fit_Cmin))
        print("\n-------------------------------------------------------------------\n")

        
        print("The covariance matrix is:\n" + str(pcov) + "\n")
        ###pg75 in theory notebook
        corr=generatecorrelation(pcov)
        print("The correlation matrix is:\n" + str(corr))
        print("You want the off-diagonal elements to be ~0, ~1 means those two fit parameters are highly coupled.\nIf they are highly coupled, you may consider adjusting your model.\n")

        #how good was the fit? Calculating Chi Squared...
        dummy=list()
        for i in range(len(y)):
            dummy.append( ( (y[i]-logisticfn(x[i],fit_T0,fit_Cmin,fit_L,fit_k)) / s[i] )**2 )
        chisq=sum(dummy)
        del dummy
        print("Chi Squared: " + str(chisq))
        chisqreduced=chisq/(len(y)-numberofparameters_logistic)
        print("Chi Squared Reduced: " + str(chisqreduced) + "\n")

        #how good was the fit? Calculating R Squared...
        ybar=sum([i for i in y])/len(y) #average of y
        SStot=sum([ (i-ybar)**2 for i in y]) #total sum of squares
        SSres=sum([ (y[i]-logisticfn(x[i],fit_T0,fit_Cmin,fit_L,fit_k))**2 for i in range(len(y))]) #sum of squares of residuals
        rsquared=1-(SSres/SStot)
        print("R Squared: " + str(rsquared))
        print(str(100*rsquared) + "% of the variation in y can be explained by the model")

        T0s.append(fit_T0)
        err_T0s.append(err_fit_T0)
        Cmins.append(fit_Cmin)
        err_Cmins.append(err_fit_Cmin)
        Ls.append(fit_L)
        err_Ls.append(err_fit_L)
        ks.append(fit_k)
        err_ks.append(err_fit_k)

        chisq_T0s.append(chisqreduced)


    return frequencies,T0s,err_T0s,Cmins,err_Cmins,Ls,err_Ls,ks,err_ks,chisq_T0s


#this function gets data for arrhenius plot
def getarrhenius(frequencies,T0s,err_T0s):
    
    #assuming no uncertainty in frequency
    arrheniusX=list();arrheniusY=list()
    err_arrheniusX=list();err_arrheniusY=list()

    for i in range(len(frequencies)):
        arrheniusX.append(1/T0s[i])
        err_arrheniusX.append( err_T0s[i] / (T0s[i]**2) )
        arrheniusY.append( log(frequencies[i] / (T0s[i]**2) ) )
        err_arrheniusY.append( 2*err_T0s[i] / T0s[i] ) #remember you have to propagate err through log

    return arrheniusX,arrheniusY,err_arrheniusX,err_arrheniusY


#this function fits arrhenius data and determines defect energy etc
def arrheniusfitting(arrheniusX,arrheniusY,err_arrheniusX,err_arrheniusY):

    print("\n\n-------------------------------------------------------------------")
    print("-------------------------------------------------------------------")
    print("Performing arrhenius fit with equation y= (b-b0) + m(x-x0)")
    print("-------------------------------------------------------------------")
    print("-------------------------------------------------------------------\n\n")

    #remember function is ln(f/T0^2) = -Ed/k * 1/T + ln(nu0)
    #                     arrheniusY = -Ed/k * arrheniusX + ln(nu0)
    initialguess=[log(1e12),-0.3/boltzmann] #guessing Ed=0.3eV and nu0=1e12
    global b0 ; global x0
    b0=sum(arrheniusY)/float(len(arrheniusY)) #average y on plot
    x0=sum(arrheniusX)/float(len(arrheniusX)) #average x on plot

    realdata = RealData(numpy.array(arrheniusX),numpy.array(arrheniusY),numpy.array(err_arrheniusX),numpy.array(err_arrheniusY)) #put our data in a form scipy.odr can use
    model = Model(linearfn) #declare our model to scipy.odr
    odr = ODR(realdata,model,initialguess) #declare which data, model, and initialguess to use in odr job
    odr.set_job(fit_type=0) #we are using explicit orthogonal distance regression (fit_type=0), we could do typical least-squares regression w/ fit_type=2, but then err in x isn't used
    output=odr.run() #do the fit
    #see all attributes of output object with dir(output)
    beta=output.beta #our parameter array
    fit_b=beta[0] #intercept part of fit
    fit_m=beta[1] #slope of fit
    
    err_fit_b=sqrt(numpy.diag(output.cov_beta)[0])
    err_fit_m=sqrt(numpy.diag(output.cov_beta)[1])
    #stderr can also be collected from output.sd_beta, but this number is slightly different. Some kind of rounding thing? Using output.cov_beta because that's the covariance matrix sort of like pcov in the curve_fit method
    pcov=output.cov_beta
    print("m = " + str(fit_m) + " +/- " + str(err_fit_m))
    print("b = " + str(fit_b) + " +/- " + str(err_fit_b))
    print("x0 = " + str(x0))
    print("b0 = " + str(b0) + "\n")

    #Now calculating parameters from the arrhenius equation.
    #remember function is ln(f/T0^2) = -Ed/k * 1/T + ln(nu0)
    print("Arrhenius fit parameters")
    Ed = -boltzmann*fit_m
    #nu0 = exp(fit_b)#this is wrong
    nu0 = exp(linearfn(beta,0))

    err_Ed = -boltzmann*err_fit_m
    err_nu0 = nu0*err_fit_b
    print("Ed = " + str(Ed) + " +/- " + str(err_Ed))
    print("nu0 = " + str('{:e}'.format(nu0)) + " +/- " + str('{:e}'.format(err_nu0))) #'{:e}'.format() converts to scientific notation

    print("") #blank line
    print("Put this equation in gnuplot etc. Using x for 1/temperature in 1/K etc")
    print("f(x) = ( " + str(fit_b) + " - " + str(b0) + " ) + " + str(fit_m) + " * ( x - " + str(x0) + " )")


    print("\n-------------------------------------------------------------------\n")

    print("The covariance matrix is:\n" + str(pcov) + "\n")
    ###pg75 in theory notebook
    corr=generatecorrelation(pcov)
    print("The correlation matrix is:\n" + str(corr))
    print("You want the off-diagonal elements to be ~0, ~1 means those two fit parameters are highly coupled.\nIf they are highly coupled, you may consider adjusting your model.\n")

    #how good was the fit? Calculating Chi Squared...
    dummy=list()
    for i in range(len(arrheniusY)):
        ################this doesn't work for uncertainty in both x and y, ask mark about this
        dummy.append( ( (arrheniusY[i]-linearfn([fit_b,fit_m],arrheniusX[i])) / err_arrheniusY[i] )**2 )
#        dummy.append( 
#        ( (arrheniusY[i]-linearfn([fit_b,fit_m],arrheniusX[i])) / err_arrheniusY[i] )**2 
#        ( (arrheniusX[i]-linearfn([fit_b,fit_m],arrheniusX[i])) / err_arrheniusX[i] )**2 
#                           this term needs to produce the x on the fit fn
#        )
    chisq=sum(dummy)
    del dummy
    print("Chi Squared: " + str(chisq))
    chisqreduced=chisq/(len(arrheniusY)-numberofparameters_linear)
    print("Chi Squared Reduced: " + str(chisqreduced) + "\n")
    chisq_arrhenius=chisqreduced

    #how good was the fit? Calculating R Squared...
    ybar=sum([i for i in arrheniusY])/len(arrheniusY) #average of y
    SStot=sum([ (i-ybar)**2 for i in arrheniusY]) #total sum of squares
    SSres=sum([ (arrheniusY[i]-linearfn([fit_b,fit_m],arrheniusX[i]))**2 for i in range(len(arrheniusY))]) #sum of squares of residuals
    rsquared=1-(SSres/SStot)
    print("R Squared: " + str(100*rsquared) + "% of the variation in y can be explained by the model")

    return Ed,nu0,err_Ed,err_nu0,chisq_arrhenius


def generatecorrelation(cov):
    corr=cov #just to mimic it's size.
    for i in range(len(cov)): #each row
        for j in range(len(cov[i])): #each column
            corr[i][j]=cov[i][j]/sqrt( numpy.diag(cov)[i] * numpy.diag(cov)[j]  )
    return corr

#this function writes data to file
def writeData(rangeNumber,Ed,nu0,err_Ed,err_nu0,chisq_arrhenius,frequencies,T0s,err_T0s,chisq_T0s,arrheniusX,arrheniusY,err_arrheniusX,err_arrheniusY):

    print("\n-------------------------------------------------------------------")
    print("-------------------------------------------------------------------\n")
    print("Outputting data to: " + locOutput)
    f=open(locOutput,'w')
    f.write("For defect " + str(rangeNumber) + ", we calculate:\n")
    f.write("Ed = " + str(Ed) + " +/- " + str(err_Ed) + " eV\n")
    f.write("nu0 = " + str(nu0) + " +/- " + str(err_nu0) + " Hz\n")
    f.write("chi squared reduced for arrhenius fit" + "\n")
    f.write("For each frequency,\n")
    f.write("T0 [K]" + delimiter + \
            "err in T0 (one stdev) [K]" + delimiter + \
            "Frequency [Hz]" + delimiter + \
            "chi squared reduced for T0 fit" + delimiter + \
            "1/T0 [1/K]" + delimiter + \
            "err in 1/T0 [1/K]" + delimiter + \
            "ln(F/T0^2) [ln(Hz/K^2)]" + delimiter + \
            "err in ln(F/T0^2) [ln(Hz/K^2)]" + delimiter + \
            "\n" \
           )
    for i in range(len(T0s)):
        f.write(str(T0s[i]) + delimiter + \
                str(err_T0s[i]) + delimiter + \
                str(frequencies[i]) + delimiter + \
                str(chisq_T0s[i]) + delimiter + \
                str(arrheniusX[i]) + delimiter + \
                str(err_arrheniusX[i]) + delimiter + \
                str(arrheniusY[i]) + delimiter + \
                str(err_arrheniusY[i]) + delimiter + \
                "\n" \
               )
    f.close()

    return



#section needed for if the user is executing this .py directly, not importing to another script or shell
if __name__ == '__main__':
    main(locInput,locOutput,stepRanges)
