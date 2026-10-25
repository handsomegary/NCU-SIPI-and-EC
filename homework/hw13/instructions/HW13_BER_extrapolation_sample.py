


#%% Copy the code of the last homework here




#%% normalize the jitter histogram to make it sum to 1 (PMF)
jtx1 = jtx1 / np.sum(jtx1)
jtx2 = 
jtx3 = 
jtx4 = 
jty1 = 
jty2 = 
jty3 = 
jty4 = 


# compute sampled jitter cdf (non-extrapolated)
cdfx1 = # hint: cumsum
cdfx2 = 
cdfx3 = 
cdfx4 = 
cdfy1 = 
cdfy2 = 
cdfy3 = 
cdfy4 = 


tj = np.linspace(-Tb/2,Tb/2, ? ) # time vector for cdf
tjp = tj # tj in ps

plt.figure(dpi=250)
plt.semilogy(tjp,cdfx1,label='cdfx1')
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.xlabel(
plt.ylabel(
plt.legend(


#%% compute sampled bathtub curve (non-extrapolated)

# First, expand the sampled cdf to 3UI wide
cdfx1_3UI = np.hstack((np.zeros(NpUI),cdfx1,np.ones(NpUI)))
cdfx2_3UI = 
cdfx3_3UI = 
cdfx4_3UI = 
cdfy1_3UI = 
cdfy2_3UI = 
cdfy3_3UI = 
cdfy4_3UI = 

# Next, compute the BER cufve
rhoT = 0.5 # assume transition density = 0.5
tBER = np.arange(2*NpUI)*dt # time vector for BER bath
tBERps = tBER*1e12

BER_x1 = rhoT*(1-cdfx1_3UI[NpUI:]) + rhoT*cdfx1_3UI[:2*NpUI]
BER_x2 = 
BER_x3 = 
BER_x4 = 
BER_y1 = 
BER_y2 = 
BER_y3 = 
BER_y4 = 

plt.figure(dpi=250)
plt.semilogy(tBERps,BER_x1,label='BERx1')
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.xlabel(
plt.ylabel(
plt.legend(



#%% convert PMF to PDF, with unit [probability/ps]
jtx1 /= 
jtx2 /= 
jtx3 /= 
jtx4 /= 
jty1 /= 
jty2 /= 
jty3 /= 
jty4 /= 



plt.figure(dpi=250)
plt.semilogy(tjp,jtx1,label='jtx1')
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.semilogy(
plt.xlabel(
plt.ylabel(
plt.legend(

#%% define function: Gaussian mixture model, PDF

# PDF
def PDF2G(t, mu1, s1, tau1, mu2, s2, tau2):
    return 

#%% define function Gaussian mixture model, CDF
from scipy.special import erfc

def CDF2G(t, 
    return 



#%% tail fitting for jty3
mu1 = 
mu2 = 
s1 =  # sigma1
s2 =  # sigma2
tau1 = 
tau2 = 


PDF2G_y3_ps = # use PDF2G function and tj as time vector


ti = Tb/2
tip1 =  # T(i+1)

BER_y3fit = # use CDF2G function and tBER, with ti and tip1

plt.figure(dpi=250,figsize=(10,4))
plt.subplot(121)
plt.semilogy(tjp,jty3,label='Sampled PDF')
plt.semilogy(tjp,PDF2G_y3_ps,'--',label='2 Gaussians')
plt.xlabel(
plt.ylabel(
plt.legend(
plt.ylim(


plt.subplot(122)
plt.semilogy(      ,label='Sampled BER'
plt.semilogy(
plt.xlabel(
plt.ylabel(
plt.grid(
plt.xlim(
plt.ylim(
plt.legend(









#%% tail fitting for jty4
# the same as for jty3
mu1 = 
mu2 = 
s1 = # sigma1
s2 = # sigma2
tau1 = 
tau2 = 


PDF2G_y4_ps = 

BER_y4fit = 




%matplotlib inline
plt.figure(dpi=




