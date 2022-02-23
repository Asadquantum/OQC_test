import numpy as np

def X(theta):
    '''X-rotation'''
    return np.array([[np.cos(theta/2),-1j*np.sin(theta/2)],[-1j*np.sin(theta/2),np.cos(theta/2)]])
def Z(theta):
    '''Z-rotation'''
    return np.array([[np.exp(-1j*theta/2),0],[0,np.exp(1j*theta/2)]])
	
def angles_of_3_Euler_pulses(U,length_Z,length_X):
    '''We now convert our unitary into three pulses with rotation angles 𝛼, 𝛽 and 𝛾'''
    x_1=np.real(U[0,0])
    x_2=np.imag(U[0,0])
    x_3=np.real(U[0,1])
    x_4=np.imag(U[0,1])

    if length_Z<length_X:
        if x_1==0:
            theta_x2_x1=np.inf
        else:
            theta_x2_x1=x_2/x_1
        if x_2==0:
            theta_x3_x2=np.inf
        else:
            theta_x3_x2=x_3/x_2
        if x_4==0:
            theta_x3_x4=np.inf
        else:
            theta_x3_x4=x_3/x_4

        alpha=np.arctan(-theta_x2_x1)+np.arctan(theta_x3_x4)
        gamma=np.arctan(-theta_x2_x1)-np.arctan(theta_x3_x4)
        beta=2*np.arctan((theta_x3_x2*np.sin((alpha+gamma)/2))/(np.sin((alpha-gamma)/2)))
        U_optimal_pulse_1=np.matmul(X(beta),Z(gamma))
        U_optimal_pulse=np.matmul(Z(alpha),U_optimal_pulse_1)
    else:
        if x_1==0:
            theta_x4_x1=np.inf
            theta_x3_x1=np.inf
        else:
            theta_x4_x1=x_4/x_1
            theta_x3_x1=x_3/x_1
        if x_2==0:
            theta_x3_x2=np.inf
        else:
            theta_x3_x2=x_3/x_2
            

        alpha=np.arctan(-theta_x4_x1)+np.arctan(-theta_x3_x2)
        gamma=np.arctan(-theta_x4_x1)-np.arctan(-theta_x3_x2)
        beta=2*np.arctan((-theta_x3_x1*np.cos((alpha+gamma)/2))/(np.sin((alpha-gamma)/2)))
        U_optimal_pulse_1=np.matmul(Z(beta),X(gamma))
        U_optimal_pulse=np.matmul(X(alpha),U_optimal_pulse_1)
    print('\n The overall rotation using optimal pulse sequence is the same: \n')
    print(U_optimal_pulse)
    return alpha,beta,gamma
	
def mainfunc(list_of_rotation_type,list_of_rotation_angle,length_Z,length_X):
    '''Main function where we first make combined unitary before dividing it into pulses'''
    ### generating combination U
    U=np.array([[1,0],[0,1]])
    for i in range(len(list_of_rotation_type)):
        if list_of_rotation_type[i]=='X':
            U_x=X(list_of_rotation_angle[i]*(np.pi/180))
            U=np.matmul(U,U_x)
        else:
            U_z1=Z(np.pi/2)
            U_x=X(list_of_rotation_angle[i]*(np.pi/180))
            U_z2=Z(-np.pi/2)
            U_y=np.matmul(U_z1,np.matmul(U_x,U_z2))
            U=np.matmul(U,U_y)
    print('\n The overall rotation using input pulse sequence is: \n')
    print(U)
    ### Identifying three pulses now
    alpha,beta,gamma= angles_of_3_Euler_pulses(U,length_Z,length_X)
    print('\n And the optimal pulse sequence is:\n')
    if length_Z<length_X: 
        print('Z ( ',alpha*(180/np.pi),' ) X (',beta*(180/np.pi),') Z (',gamma*(180/np.pi),')')
    else: 
        print('X ( ',alpha*(180/np.pi),' ) Z (',beta*(180/np.pi),') X (',gamma*(180/np.pi),')')
    