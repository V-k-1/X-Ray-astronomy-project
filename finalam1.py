import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import sqrt
def master(x_ar,a1=0,c1=0,s1=1,a2=0,c2=0,s2=1,a3=0,c3=0,s3=1,a4=0,c4=0,s4=1,a5=0,c5=0,s5=1,a6=0,c6=0,s6=1,a7=0,c7=0,s7=1,a8=0,c8=0,s8=1,a9=0,c9=0,s9=1,a10=0,c10=0,s10=1,a11=0,c11=0,s11=1):
    return [_1gaussian(i,a1,c1,s1)+_1gaussian(i,a2,c2,s2)+_1gaussian(i,a3,c3,s3)+_1gaussian(i,a4,c4,s4)+_1gaussian(i,a5,c5,s5)+_1gaussian(i,a6,c6,s6)+_1gaussian(i,a7,c7,s7)+_1gaussian(i,a8,c8,s8)+_1gaussian(i,a9,c9,s9)+_1gaussian(i,a10,c10,s10)+_1gaussian(i,a11,c11,s11) for i in x_ar]
def _1gaussian(x, a1, c1, s1):
    return (a1*(1/((s1)*(np.sqrt(2*np.pi))))*(np.exp(-((x-c1)**2/((2*(s1))**2)))))
def fwhm(x):
    return 2*sqrt(2*np.log(2))*(x)
def straightline(x,m,c):
    return x*m+c
def energy(x):
    return calfac*(x)+inter
def plot1g(x, var):
    for i in range(int(len(var)/3)):
        print(*list(var[i*3+j] for j in range(3)))
        plt.plot(x, _1gaussian(x, *list(var[i*3+j] for j in range(3))),'--',lw=0.6)
a=np.loadtxt('am1.txt')
b=np.loadtxt('background.txt')
a[:,1]=a[:,1]-b[:,1]
for i in range(len(a[:,1])):
    if(i<90) : a[i,1]=0

x_ar,y_ar,maxch,maxc,I=a[:1000,0],a[:1000,1],[],[],[]
plt.plot(x_ar,y_ar,'g')

for i in range(len(x_ar)):
    if(y_ar[i]==max(y_ar[i-min(1,i):i+1])):
        if(y_ar[i]>10) and (y_ar[i]-y_ar[i-3])>4 and y_ar[i]-y_ar[i+3]>4:
            maxc.append(y_ar[i])
            maxch.append(x_ar[i])
print("max. count values-",maxc)
print("max. count at channel values-",maxch)

calfac=(59.54-13.94)/(928.4730782104126-217.5)
inter=59.54-calfac*928.4730782104126
ee=-1*inter/calfac
print("calibration factor-",calfac,"intercept-",inter)
E=[]
for i in range(len(maxch)):
    E.append(calfac*maxch[i])
print("energy values according to peaks-",E)

x1=[94,181,211,277,320,406,492,560,925]
x2=[109,195,228,286,341,418,527,576,934]
#C=[[100,102,104],[183,186,188,190,193],[218,220,223,0],[278,280,0],[325,327,335,0,0],[407,409,413,416],[562,566,568,571,0],[926,928,0]]
C=[[0],[0],[0],[0],[0],[0],[502,519],[0],[0]]
A=[]
channel,E,ER,FWHM,channelerr,FWHMerr,ERerr,Eerr=[],[],[],[],[],[],[],[]
for i in range(len(x1)):
    print("peak ",i+1)
    x_ar=a[int(x1[i]):int(x2[i]),0]
    y_ar=a[int(x1[i]):int(x2[i]),1]
    p0=[]
    bound=[],[]
    for k in range(len(C[i])):
        p0.append(max(y_ar))
        if C[i][k]!=0 :
           p0.append(C[i][k])
        else :
           p0.append(sum(x_ar)/len(x_ar))
        p0.append(0.2)
        bound[0].append(0)
        bound[1].append(10000000)
        if C[i][k]!=0 :
            bound[0].append(C[i][k]-0.5)
            bound[1].append(C[i][k]+0.5)
        else :
            bound[0].append(x_ar[0])
            bound[1].append(x_ar[-1])
        bound[0].append(0)
        bound[1].append(10)
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
        Eerr.append(energy(ee+pcov_2gauss[3*j+1,3*j+1]**0.5))
        FWHM.append(fwhm(popt_2gauss[3*j+2]))
        FWHMerr.append(fwhm(pcov_2gauss[3*j+2,3*j+2]**0.5))
        ER.append((fwhm(popt_2gauss[3*j+2])/(popt_2gauss[3*j+1]))*100)
        ERerr.append(ER[-1]*(sqrt((FWHMerr[-1]/FWHM[-1])**2+(Eerr[-1]/E[-1])**2)))
        print("FWHM",j+1,"-",FWHM[-1],"(+/-)",FWHMerr[-1],"\nEnergy",j+1,"-",E[-1],"(+/-)",Eerr[-1],"\nEnergy resolution",j+1,"-",ER[-1],"(+/-)",ERerr[-1],'%')
plt.show()
