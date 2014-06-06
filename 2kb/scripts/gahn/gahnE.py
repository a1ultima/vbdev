H=0.01
h=0.99
e=0.7
T0=0.0
T1=0.0
T2=0.0
T3=1.0

in1 = [H,h,T0,T1,T2,T3,e]


#1*H*H*T1*T1+e*2*H*h*T1*T1
#2*H*H*T1*T2+e*4*H*h*T1*T2
#2*H*H*T1*T3+e*4*H*h*T1*T3
#1*H*H*T2*T2+e*2*H*h*T2*T2
#2*H*H*T2*T3+e*4*H*h*T2*T3
#1*H*H*T3*T3+e*2*H*h*T3*T3
#2*(1-e)*H*h*T1*T1
#4*(1-e)*H*h*T1*T2
#4*(1-e)*H*h*T1*T3
#2*(1-e)*H*h*T2*T2
#4*(1-e)*H*h*T2*T3
#2*(1-e)*H*h*T3*T3
#1*h*h*T1*T1
#2*h*h*T1*T2
#2*h*h*T1*T3
#1*h*h*T2*T2
#2*h*h*T2*T3
#1*h*h*T3*T3

def out( in1 ):

	H		=in1[0]
	h		=in1[1]
	T0		=in1[2]
	T1		=in1[3]
	T2		=in1[4]
	T3		=in1[5]
	e		=in1[6]

	TotalF = (1*H*H*T1*T1+e*2*H*h*T1*T1+2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T1*T3+e*4*H*h*T1*T3+1*H*H*T2*T2+e*2*H*h*T2*T2+2*H*H*T2*T3+e*4*H*h*T2*T3+1*H*H*T3*T3+e*2*H*h*T3*T3+2*(1-e)*H*h*T1*T1+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+2*(1-e)*H*h*T2*T2+4*(1-e)*H*h*T2*T3+2*(1-e)*H*h*T3*T3+1*h*h*T1*T1+2*h*h*T1*T2+2*h*h*T1*T3+1*h*h*T2*T2+2*h*h*T2*T3+1*h*h*T3*T3)
	#H_2 = (1*(1*H*H*T1*T1+e*2*H*h*T1*T1 + 2*H*H*T1*T2+e*4*H*h*T1*T2 + 2*H*H*T1*T3+e*4*H*h*T1*T3 + 1*H*H*T2*T2+e*2*H*h*T2*T2 + 2*H*H*T2*T3+e*4*H*h*T2*T3 + 1*H*H*T3*T3+e*2*H*h*T3*T3) + 0.5*(2*(1-e)*H*h*T1*T1 + 4*(1-e)*H*h*T1*T2 + 4*(1-e)*H*h*T1*T3 + 2*(1-e)*H*h*T2*T2 + 4*(1-e)*H*h*T2*T3 + 2*(1-e)*H*h*T3*T3))/TotalF
	H_2 = H*(1+e*h)
	h_2 = 1-H_2 # (1*(1*h*h*T1*T1+2*h*h*T1*T2+2*h*h*T1*T3+1*h*h*T2*T2+2*h*h*T2*T3+1*h*h*T3*T3)+0.5*(2*(1-e)*H*h*T1*T1+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+2*(1-e)*H*h*T2*T2+4*(1-e)*H*h*T2*T3+2*(1-e)*H*h*T3*T3))/TotalF
	T0_2= (1*(1*H*H*T1*T1+e*2*H*h*T1*T1) + 0.5*(2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T1*T3+e*4*H*h*T1*T3))/TotalF
	T1_2= (1*(1*H*H*T2*T2+e*2*H*h*T2*T2+2*(1-e)*H*h*T1*T1+1*h*h*T1*T1)+0.5*(2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T2*T3+e*4*H*h*T2*T3+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+2*h*h*T1*T2+2*h*h*T1*T3))/TotalF
	T2_2= (1*(1*H*H*T3*T3+e*2*H*h*T3*T3+2*(1-e)*H*h*T2*T2+1*h*h*T2*T2)+0.5*(2*H*H*T1*T3+e*4*H*h*T1*T3+2*H*H*T2*T3+e*4*H*h*T2*T3+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T2*T3+2*h*h*T1*T2+2*h*h*T2*T3))/TotalF
	T3_2= 1-(T0_2+T1_2+T2_2)
	return [H_2, h_2, T0_2, T1_2, T2_2, T3_2, e]

print('gahn is cool: ')

generationList= []
for i in range(0,10):
	print('\tGeneration: '+str(i))
	in1 = out(in1)
	generationList.append(in1)

import numpy as np 

moo = np.array(generationList).T
 
import pprint
pprint.pprint(moo[0])
