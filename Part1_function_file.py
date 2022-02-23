import numpy as np
def X(theta):
    '''X-rotation'''
    return np.array([[np.cos(theta/2),-1j*np.sin(theta/2)],[-1j*np.sin(theta/2),np.cos(theta/2)]])
def Y(theta):
    '''Y-rotation'''
    return np.array([[np.cos(theta/2),-np.sin(theta/2)],[np.sin(theta/2),np.cos(theta/2)]])
	
def angles_of_3_Euler_pulses(U):
    '''We now convert our unitary into three pulses with rotation angles 𝛼, 𝛽 and 𝛾'''
    x_1=np.real(U[0,0])
    x_2=np.imag(U[0,0])
    x_3=np.real(U[0,1])
    x_4=np.imag(U[0,1])

    if x_1==0:
        theta_x4_x1=np.inf
        theta_x3_x1=np.inf
    else:
        theta_x4_x1=x_4/x_1
        theta_x3_x1=x_3/x_1
    if x_3==0:
        theta_x2_x3=np.inf
    else:
        theta_x2_x3=x_2/x_3

    alpha=np.arctan(-theta_x4_x1)+np.arctan(theta_x2_x3)
    gamma=np.arctan(-theta_x4_x1)-np.arctan(theta_x2_x3)
    beta=2*np.arctan((-theta_x3_x1*np.cos((alpha+gamma)/2))/(np.cos((alpha-gamma)/2)))
    U_optimal_pulse_1=np.matmul(Y(beta),X(gamma))
    U_optimal_pulse=np.matmul(X(alpha),U_optimal_pulse_1)
    print('\n The overall rotation using optimal pulse sequence is the same: \n')
    print(U_optimal_pulse)
    return alpha,beta,gamma
	
def mainfunc(list_of_rotation_type,list_of_rotation_angle):
    '''Main function where we first make combined unitary before dividing it into pulses'''
    ### generating combination U
    U=np.array([[1,0],[0,1]])
    for i in range(len(list_of_rotation_type)):
        if list_of_rotation_type[i]=='X':
            U_x=X(list_of_rotation_angle[i]*(np.pi/180))
            U=np.matmul(U,U_x)
        else:
            U_y=Y(list_of_rotation_angle[i]*(np.pi/180))
            U=np.matmul(U,U_y)
    print('\n The overall rotation using input pulse sequence is: \n')
    print(U)
    ### Identifying three pulses now
    alpha,beta,gamma= angles_of_3_Euler_pulses(U)
    print('\n And the optimal pulse sequence is:\n')
    print('X ( ',alpha*(180/np.pi),' ) Y (',beta*(180/np.pi),') X (',gamma*(180/np.pi),')')