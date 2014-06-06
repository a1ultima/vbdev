
	H		=in1[0]
	h		=in1[1]
	T0		=in1[2]
	T1		=in1[3]
	T2		=in1[4]
	T3		=in1[5]
	T4		=in1[6]

1*H*H*T1*T1+e*2*H*h*T1*T1
2*H*H*T1*T2+e*4*H*h*T1*T2
2*H*H*T1*T3+e*4*H*h*T1*T3
2*H*H*T1*T4+e*4*H*h*T1*T4
1*H*H*T2*T2+e*2*H*h*T2*T2
2*H*H*T2*T3+e*4*H*h*T2*T3
2*H*H*T2*T4+e*4*H*h*T2*T4
1*H*H*T3*T3+e*2*H*h*T3*T3
2*H*H*T3*T4+e*4*H*h*T3*T4
2*(1-e)*H*h*T1*T1
4*(1-e)*H*h*T1*T2
4*(1-e)*H*h*T1*T3
4*(1-e)*H*h*T1*T4
2*(1-e)*H*h*T2*T2
4*(1-e)*H*h*T2*T3
4*(1-e)*H*h*T2*T4
2*(1-e)*H*h*T3*T3
4*(1-e)*H*h*T3*T4
1*h*h*T1*T1
2*h*h*T1*T2
2*h*h*T1*T3
2*h*h*T1*T4
1*h*h*T2*T2
2*h*h*T2*T3
2*h*h*T2*T4
1*h*h*T3*T3
1*h*h*T3*T4

TotalSurvivors 	= (1*H*H*T1*T1+e*2*H*h*T1*T1+2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T1*T3+e*4*H*h*T1*T3+2*H*H*T1*T4+e*4*H*h*T1*T4+1*H*H*T2*T2+e*2*H*h*T2*T2+2*H*H*T2*T3+e*4*H*h*T2*T3+2*H*H*T2*T4+e*4*H*h*T2*T4+1*H*H*T3*T3+e*2*H*h*T3*T3+2*H*H*T3*T4+e*4*H*h*T3*T4+1*H*H*T4*T4+e*2*H*h*T4*T4+2*(1-e)*H*h*T1*T1+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+4*(1-e)*H*h*T1*T4+2*(1-e)*H*h*T2*T2+4*(1-e)*H*h*T2*T3+4*(1-e)*H*h*T2*T4+2*(1-e)*H*h*T3*T3+4*(1-e)*H*h*T3*T4+2*(1-e)*H*h*T4*T4+1*h*h*T1*T1+2*h*h*T1*T2+2*h*h*T1*T3+2*h*h*T1*T4+1*h*h*T2*T2+2*h*h*T2*T3+2*h*h*T2*T4+1*h*h*T3*T3+2*h*h*T3*T4+1*h*h*T4*T4)
TotalDead 		= 1 - TotalSurvivors
Hnext			= (1*(1*H*H*T1*T1+e*2*H*h*T1*T1+2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T1*T3+e*4*H*h*T1*T3+2*H*H*T1*T4+e*4*H*h*T1*T4+1*H*H*T2*T2+e*2*H*h*T2*T2+2*H*H*T2*T3+e*4*H*h*T2*T3+2*H*H*T2*T4+e*4*H*h*T2*T4+1*H*H*T3*T3+e*2*H*h*T3*T3+2*H*H*T3*T4+e*4*H*h*T3*T4+1*H*H*T4*T4+e*2*H*h*T4*T4)+(1/2)*(2*(1-e)*H*h*T1*T1+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+4*(1-e)*H*h*T1*T4+2*(1-e)*H*h*T2*T2+4*(1-e)*H*h*T2*T3+4*(1-e)*H*h*T2*T4+2*(1-e)*H*h*T3*T3+4*(1-e)*H*h*T3*T4+2*(1-e)*H*h*T4*T4))/TotalSurvivors
hnext 			= (1*(1*h*h*T1*T1+2*h*h*T1*T2+2*h*h*T1*T3+2*h*h*T1*T4+1*h*h*T2*T2+2*h*h*T2*T3+2*h*h*T2*T4+1*h*h*T3*T3+2*h*h*T3*T4+1*h*h*T4*T4)+(1/2)*(2*(1-e)*H*h*T1*T1+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+4*(1-e)*H*h*T1*T4+2*(1-e)*H*h*T2*T2+4*(1-e)*H*h*T2*T3+4*(1-e)*H*h*T2*T4+2*(1-e)*H*h*T3*T3+4*(1-e)*H*h*T3*T4+2*(1-e)*H*h*T4*T4))/TotalSurvivors
hnext 			= 1-Hnext
T0next 			= (1*(1*H*H*T1*T1+e*2*H*h*T1*T1)+(1/2)*(2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T1*T3+e*4*H*h*T1*T3+2*H*H*T1*T4+e*4*H*h*T1*T4))/TotalSurvivors
T1next			= (1*(1*H*H*T2*T2+e*2*H*h*T2*T2+2*(1-e)*H*h*T1*T1+1*h*h*T1*T1)+(1/2)*(2*H*H*T1*T2+e*4*H*h*T1*T2+2*H*H*T2*T3+e*4*H*h*T2*T3+2*H*H*T2*T4+e*4*H*h*T2*T4+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T1*T3+4*(1-e)*H*h*T1*T4+2*h*h*T1*T2+2*h*h*T1*T3+2*h*h*T1*T4))/TotalSurvivors
T2next			= (1*(1*H*H*T3*T3+e*2*H*h*T3*T3+2*(1-e)*H*h*T2*T2+1*h*h*T2*T2)+(1/2)*(2*H*H*T1*T3+e*4*H*h*T1*T3+2*H*H*T2*T3+e*4*H*h*T2*T3+2*H*H*T3*T4+e*4*H*h*T3*T4+4*(1-e)*H*h*T1*T2+4*(1-e)*H*h*T2*T3+4*(1-e)*H*h*T2*T4+2*h*h*T1*T2+2*h*h*T2*T3+2*h*h*T2*T4))/TotalSurvivors
T3next			= (1*(1*H*H*T4*T4+e*2*H*h*T4*T4+2*(1-e)*H*h*T3*T3+1*h*h*T3*T3)+(1/2)*(2*H*H*T1*T4+e*4*H*h*T1*T4+2*H*H*T2*T4+e*4*H*h*T2*T4+2*H*H*T3*T4+e*4*H*h*T3*T4+4*(1-e)*H*h*T1*T3+4*(1-e)*H*h*T2*T3+4*(1-e)*H*h*T3*T4+2*h*h*T1*T3+2*h*h*T2*T3+2*h*h*T3*T4))/TotalSurvivors
T4next			= (1*(2*(1-e)*H*h*T4*T4+1*h*h*T4*T4)+(1/2)*(4*(1-e)*H*h*T1*T4+4*(1-e)*H*h*T2*T4+4*(1-e)*H*h*T3*T4+2*h*h*T1*T4+2*h*h*T2*T4+2*h*h*T3*T4))/TotalSurvivors
