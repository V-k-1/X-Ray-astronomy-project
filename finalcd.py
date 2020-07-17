#This is the program used for fitting the Cadmium spectra,the other two programs for Iron & Americium are almost similar with changes indicated 
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import sqrt
'''   declaration of functions : '_1gaussian' is funtion for a gaussian curve(parameters-amplitude, mean, std. deviation
  'master' is the function for multiple gaussian curves added together with a background noise 'd'
   plot1g helps in plotting the indivisual gaussian curves from the fitting parameters obtained from 'master'   '''
def master(x_ar,d=0,a1=0,c1=0,s1=1,a2=0,c2=0,s2=1,a3=0,c3=0,s3=1,a4=0,c4=0,s4=1,a5=0,c5=0,s5=1,a6=0,c6=0,s6=1,a7=0,c7=0,s7=1,a8=0,c8=0,s8=1,a9=0,c9=0,s9=1,a10=0,c10=0,s10=1,a11=0,c11=0,s11=1):
    return [_1gaussian(i,a1,c1,s1)+_1gaussian(i,a2,c2,s2)+_1gaussian(i,a3,c3,s3)+_1gaussian(i,a4,c4,s4)+_1gaussian(i,a5,c5,s5)+_1gaussian(i,a6,c6,s6)+_1gaussian(i,a7,c7,s7)+_1gaussian(i,a8,c8,s8)+_1gaussian(i,a9,c9,s9)+_1gaussian(i,a10,c10,s10)+_1gaussian(i,a11,c11,s11)+d for i in x_ar]
def _1gaussian(x, a1, c1, s1):
    return (a1*(1/((s1)*(np.sqrt(2*np.pi))))*(np.exp(-((x-c1)**2/((2*(s1))**2)))))
def energy(x):
    return calfac*(x)+intercept
def fwhm(x):
    return 2*sqrt(2*np.log(2))*(x)
def straightline(x,m,c):
    return [x[i]*m+c for i in range(len(x))]
def plot1g(x, var):
    for i in range(int(len(var)/3)):
        r=list([*list(var[i*3+j+1] for j in range(3)),var[0]]); print(r)
#        f.write(str(str(r)+'\n'))
        plt.plot(x, var[0]+_1gaussian(x, *list(var[i*3+j+1] for j in range(3))),'--',lw=0.6)
a=np.loadtxt('cumu_CdCu.txt')
b=np.loadtxt('cumu_blank.txt')
a[:,1]=a[:,1]-b[:,1]
for i in range(1):
    a[i,1]=0   #there's too much noise till channel 90; hence those values were deleted
plt.plot(a[:,0],a[:,1],'g')

x_ar,y_ar=a[:,0],a[:,1]
#two peaks were identified from the plots with a known definite energy, and their means were used to compute the calibration factor
calfac=(24.9427-22.16317)/(437.99999999999994-388.5043402822684)
intercept=24.9427-calfac*437.99999999999994 # if intercept is not equal to zero then it means 0th channel does not correspond to 0keV energy
print("calibration factor-",calfac,"intercept-",intercept)
#specify the channel region which you want to fit(from x1[i] to x2[i])
#x1=[131,378,427,1000,1124,1154,1534]
#x2=[183,398,450,1052,1148,1187,1550]
#C[i] array should contain the(exact) channel(s) where you think there must be a peak, use 0 to specify an unknown peak,len(C[i]) will be the number of gaussian functions the program is going to fit the region with
#for eg. use C=[[0],[0],[0],[0],[0],[0],[0],[0]] to fit every region with 1 gaussian function without specifying where the peak might be
#C=[[147,175],[389],[437],[1013,1038],[1135],[1174],[0]]
x1=[42,368]
x2=[125,390]
C=[[53,93],[0]]
channel,E,ER,FWHM,channelerr,FWHMerr,ERerr,Eerr=[],[],[],[],[],[],[],[]
for i in range(len(x1)):
    print("peak ",i+1)
    x_ar=a[int(x1[i]):int(x2[i]),0]
    y_ar=a[int(x1[i]):int(x2[i]),1]
    p0=[0]
    bound=[0],[5]#[min(a[int(x1[i])-20:int(x2[i])+20,1])+2]
    for k in range(len(C[i])):
        p0.append(max(y_ar))
        if C[i][k]!=0 : p0.append(C[i][k])
        else : p0.append(sum(x_ar)/len(x_ar))
        p0.append(0.2)
        bound[0].append(0)
        bound[1].append(1000000)
        if C[i][k]!=0 :
            bound[0].append(C[i][k]-2)
            bound[1].append(C[i][k]+2)
        else :
            bound[0].append(x_ar[0])
            bound[1].append(x_ar[-1])
        bound[0].append(0)
        bound[1].append(50)
    print(p0,bound)
    #fitvars contains the fitting parameters; fiterr contains the covariance matrix
    fitvars,fiterr=scipy.optimize.curve_fit(master,x_ar,y_ar,p0,bounds=bound)
    plot1g(x_ar,fitvars)
    plt.plot(x_ar,master(x_ar,*fitvars))
    for j in range(len(C[i])):
        E.append(energy(fitvars[3*j+2]))
        channel.append(fitvars[3*j+2])
        channelerr.append(fiterr[3*j+2,3*j+2]**0.5)
        Eerr.append(energy(fiterr[3*j+2,3*j+2]**0.5))
        FWHM.append(fwhm(fitvars[3*j+3]))
        FWHMerr.append(fwhm(fiterr[3*j+3,3*j+3]**0.5))
        ER.append((fwhm(fitvars[3*j+3])/(fitvars[3*j+2]))*100)
        ERerr.append(ER[-1]*(sqrt((FWHMerr[-1]/FWHM[-1])**2+(Eerr[-1]/E[-1])**2)))
        print("FWHM",j+1,"-",FWHM[-1],"(+/-)",FWHMerr[-1],"\nEnergy",j+1,"-",E[-1],"(+/-)",Eerr[-1],"\nEnergy resolution",j+1,"-",ER[-1],"(+/-)",ERerr[-1],'%')
plt.show()

