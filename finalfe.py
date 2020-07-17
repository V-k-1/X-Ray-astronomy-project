import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import sqrt
def master(x_ar,a1=0,c1=0,s1=1,a2=0,c2=0,s2=1,a3=0,c3=0,s3=1,a4=0,c4=0,s4=1,a5=0,c5=0,s5=1,a6=0,c6=0,s6=1,a7=0,c7=0,s7=1,a8=0,c8=0,s8=1,a9=0,c9=0,s9=1,a10=0,c10=0,s10=1,a11=0,c11=0,s11=1):
    return [_1gaussian(i,a1,c1,s1)+_1gaussian(i,a2,c2,s2)+_1gaussian(i,a3,c3,s3)+_1gaussian(i,a4,c4,s4)+_1gaussian(i,a5,c5,s5)+_1gaussian(i,a6,c6,s6)+_1gaussian(i,a7,c7,s7)+_1gaussian(i,a8,c8,s8)+_1gaussian(i,a9,c9,s9)+_1gaussian(i,a10,c10,s10)+_1gaussian(i,a11,c11,s11) for i in x_ar]
def _1gaussian(x, a1, c1, s1):
    return (a1*(1/((s1)*(np.sqrt(2*np.pi))))*(np.exp(-((x-c1)**2/((2*(s1))**2)))))
def energy(x):
#    print("energy",calfac*(x-1535)+88.04)
    return calfac*(x)+intercept
def fwhm(x):
    return 2*sqrt(2*np.log(2))*(x)
def straightline(x,m,c):
    return x*m+c
def plot1g(x, var):
    for i in range(int(len(var)/3)):
        print(*list(var[i*3+j] for j in range(3)))
        plt.plot(x, _1gaussian(x, *list(var[i*3+j] for j in range(3))),'--',lw=0.6)
a=np.loadtxt('fe1.txt')
b=np.loadtxt('background.txt')
a[:,1]=a[:,1]-b[:,1]
for i in range(len(a[:,1])):
    if(i<90) : a[i,1]=0
plt.plot(a[:1000,0],a[:1000,1],'g')
x_ar,y_ar,maxch,maxc,I=a[:,0],a[:,1],[],[],[]
calfac=(5.89875-6.49045)/(106.30739009815892-116.92053364468042)
intercept=5.89875-calfac*106.30739009815892
print("calibration factor-",calfac,'intercept-',intercept)
x1=[100,199]
x2=[128,222]
channel,E,ER,FWHM,channelerr,FWHMerr,ERerr,Eerr=[],[],[],[],[],[],[],[]
for i in range(2):
    print("peak ",i+1)
    x_ar=a[int(x1[i]):int(x2[i]),0]
    y_ar=a[int(x1[i]):int(x2[i]),1]
    if i==0:
        p0=[640000,106,3,1000,117,3]
        bound=[300000,105,2,0,115,2],[900000,107,5,100000,119,5]
    else:
        p0=[1000,206,0.2,1000,208,0.2,0,203,0.2,0,201,0.2,0,210,0.2,0,211,0.2,0,213,0.2]
        bound=[0,205,0,0,207,0,0,202,0,0,200.5,0,0,209,0,0,210,0.2,0,212,0.2],[1000,207,1,1000,209,1,1000,204,1,1000,201.5,1,1000,211,1,1000,212,1,1000,214,1]
    print(p0,bound)
    popt_2gauss,pcov_2gauss=scipy.optimize.curve_fit(master,x_ar,y_ar,p0,bounds=bound)
    plot1g(x_ar,popt_2gauss)
    plt.plot(x_ar,master(x_ar,*popt_2gauss))
    b=int((len(popt_2gauss))/3)
    maxchannel=0
    maxchannelerr=0
    for j in range(b):
        E.append(energy(popt_2gauss[3*j+1]))
        channel.append(popt_2gauss[3*j+1])
        channelerr.append(pcov_2gauss[3*j+1,3*j+1]**0.5)
        Eerr.append(energy(pcov_2gauss[3*j+1,3*j+1]**0.5))
        FWHM.append(fwhm(popt_2gauss[3*j+2]))
        FWHMerr.append(fwhm(pcov_2gauss[3*j+2,3*j+2]**0.5))
        ER.append((fwhm(popt_2gauss[3*j+2])/(popt_2gauss[3*j+1]))*100)
        ERerr.append(ER[-1]*(sqrt((FWHMerr[-1]/FWHM[-1])**2+(Eerr[-1]/E[-1])**2)))
        print("FWHM",j+1,"-",FWHM[-1],"(+/-)",FWHMerr[-1],"\nEnergy",j+1,"-",E[-1],"(+/-)",Eerr[-1],"\nEnergy resolution",j+1,"-",ER[-1],"(+/-)",ERerr[-1],'%')
plt.show()
