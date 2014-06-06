import numpy as np 

# H  = .01
# h  = .99
# T0 = .0
# T1 = .0
# T2 = .0
# T3 = 1.0
# e  = 0.7
# H_new   = (H+e*h*H);
# T0_new  = ((H*(2*e*h+H)*T1)/(T1+T2+T3));
# T1_new  = ((h*(h-4*e*H)*T1+H*(2*e*h+H)*T2)/(T1+T2+T3))
# T2_new  = ((h*(h-4*e*H)*T2+H*(2*e*h+H)*T3)/(T1+T2+T3))
# T3_new  = ((h*T3*(h*(T1+T2+T3)+2*H*(T2+T3-e*(2*T1+T2+T3))))/(T1+T2+T3)^2);

def H_next( H, h, T0, T1, T2, T3, e):
    H_new   = (H+e*h*H)
    h_new   = (1.-H_new)
    T0_new  = ((H*(2.*e*h+H)*T1)/(T1+T2+T3))
    T1_new  = ((h*(h-4.*e*H)*T1+H*(2.*e*h+H)*T2)/(T1+T2+T3))
    T2_new  = ((h*(h-4.*e*H)*T2+H*(2.*e*h+H)*T3)/(T1+T2+T3))
    T3_new  = ((h*T3*(h*(T1+T2+T3)+2.*H*(T2+T3-e*(2.*T1+T2+T3))))/np.power((T1+T2+T3),2.))
    return H_new, h_new, T0_new, T1_new, T2_new, T3_new
#H, h, T0, T1, T2, T3 = H_next(.01,.99,.0,.0,.0,1.,0.7)

H   = .01
h   = .99
T0  = .0
T1  = .0
T2  = .0
T3  = 1.0
e   = 0.1
t   = 100

t_timeSeries= range(0,t)
H_timeSeries= [[H, h, T0, T1, T2, T3]]
for i in t_timeSeries[:-1]:
    H, h, T0, T1, T2, T3 = H_next(H, h, T0, T1, T2, T3, e)
    H_timeSeries.append([H, h, T0, T1, T2, T3])
H_timeSeries = np.clip(np.array(H_timeSeries).T,0,1)

# H_t  = H_timeSeries[0] 
# h_t  = H_timeSeries[1] 
# T0_t = H_timeSeries[2] 
# T1_t = H_timeSeries[3] 
# T2_t = H_timeSeries[4] 
# T3_t = H_timeSeries[5] 

import matplotlib.pyplot as plt 
# 255-255-0
# 34-139-34
# 255-140-0
# 25-25-112
# 0-0-255
# 65-105-225
# 135-206-250
# 0-206-209
# 176-196-222
# 100-149-237
# 0-191-255
# 0-255-255
# colors = [
#         ((176./255.),(196./255.),(222./255.)),  # T0 light
#         #((135./255.),(206./255.),(250./255.)), # 
#         ((0./255.),(191./255.),(255./255.)),    # T1
#         #((100./255.),(149./255.),(237./255.)), # 
#         #((0./255.),(206./255.),(209./255.)),   # 
#         #((65./255.),(105./255.),(225./255.)),  # 
#         #((0./255.),(191./255.),(255./255.)),   # 
#         ((0./255.),(0./255.),(255./255.)),      # T2
#         ((25./255.),(25./255.),(112./255.))     # T3 dark
#         ]
colors = [
        ((255./255.),(0./255.),(0./255.)),      # H  red
        ((34./255.),(139./255.),(34./255.)),    # h  green
        ((176./255.),(196./255.),(222./255.)),  # T0 light
        #((135./255.),(206./255.),(250./255.)), # 
        ((0./255.),(191./255.),(255./255.)),    # T1
        #((100./255.),(149./255.),(237./255.)), # 
        #((0./255.),(206./255.),(209./255.)),   # 
        #((65./255.),(105./255.),(225./255.)),  # 
        #((0./255.),(191./255.),(255./255.)),   # 
        ((0./255.),(0./255.),(255./255.)),      # T2
        ((25./255.),(25./255.),(112./255.))     # T3 dark
        ]
# colors = [  
#         ((255./255.),(0./255.),(0./255.)),      # red
#         ((255./255.),(140./255.),(0./255.)),    # orange
#         ((255./255.),(255./255.),(0./255.)),    # yellow
#         #((165./255.),(42./255.),(42./255.)),    # brown
#         #((0./255.),(0./255.),(0./255.))         # black
#         ((34./255.),(139./255.),(34./255.))         # green
#         ]

# dynamics = H_timeSeries[2:]
# for count,allele in enumerate(dynamics):
#     print(colors[count])
#     plt.plot(t_timeSeries,allele,color=colors[count],linewidth=2.)
# plt.ylabel('Frequency')
# plt.xlabel('Generation')
# plt.title('Telomere Length Dynamics')
# plt.show()

dynamics = H_timeSeries[:]
for count,allele in enumerate(dynamics):
    print(colors[count])
    plt.plot(t_timeSeries,allele,color=colors[count],linewidth=2.)
plt.ylabel('Frequency')
plt.xlabel('Generation')
plt.title('Telomere Length Dynamics')
plt.show()

