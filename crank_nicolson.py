import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------
#                                                       INPUT
# ----------------------------------------------------------------------------------------------------------------------

k=float(input('Enter diffusivity(k)(in cm^2/s): '))
x=float(input('Enter length of rod(in cm): '))
dx=float(input('Enter length step size: '))
t=float(input('Enter simulation time(in s): '))
dt=float(input('Enter time step size: '))

#Calculating number of nodes based on step size and length(or time)

length_nodes=int(x/dx)-1
time_nodes=int(t/dt)
alpha= k*dt/(dx**2)
# print(alpha)

#calculating number of rows and columns for the grid which will show temperature at a node n at time t

m=time_nodes+1
n=length_nodes+2
grid=np.zeros((m,n))

#initialising the grid with initial and boundary conditions

print('Enter boundary conditions(in celcius)')
t1=int(input('Enter surface 1 temperature: '))
t2=int(input('Enter surface 2 temperature: '))

for i in range(m):
    grid[i,0]=t1
    grid[i,length_nodes+1]=t2
t_initial=int(input('Enter initial rod temperature: '))
for i in range(1,length_nodes+1):
    grid[time_nodes,i]=t_initial
#print(grid)

#Forming the tridiagonal array based on alpha value

diag=2*(1+alpha)
below_diag=-alpha
above_diag=-alpha

left_matrix=np.zeros((length_nodes,length_nodes))
for i in range(length_nodes-1):
            left_matrix[i][i]=diag
            left_matrix[i][i+1]=above_diag
            left_matrix[i+1][i]=below_diag
left_matrix[length_nodes-1][length_nodes-1]=diag


# -------------------------------------------------------------------------------------------------------------------------
#                                                      CALCULATION
# ------------------------------------------------------------------------------------------------------------------------





temp=np.zeros((length_nodes,1))                         #temperature matrix(unknown)
coeff_matrix=np.zeros((length_nodes,1))                 #coefficient matrix calculated using known temperatures on grid

#based on the position of node the coefficient matrix is calculated

for i in range(time_nodes,0,-1):
    for j in range(1,length_nodes+1,1):
        #left node
        if j==1:
            U=grid[i,j]                                     #temperature at current node
            Uw=grid[i,j-1]                                  #temperature at node to the west
            Ue=grid[i,j+1]                                  #temperature at node to the east
            Unw=grid[i-1,j-1]                               #temperature at node to the north-west
            Une=grid[i-1,j+1]                               #temperature at node to the north-east
            coeff_matrix[j-1,0]=(alpha*(Uw+Unw+Ue))+(2*(1-alpha)*U)    
           
       #right node
        elif j==(length_nodes):
            U=grid[i,j]
            Uw=grid[i,j-1]
            Ue=grid[i,j+1]
            Une=grid[i-1,j+1]
            coeff_matrix[j-1,0]=(alpha*(Uw+Ue+Une))+(2*(1-alpha)*U)
           
         
        #center node
        else:
            U=grid[i,j]
            Uw=grid[i,j-1]
            Ue=grid[i,j+1]
            coeff_matrix[j-1,0]=(alpha*(Uw+Ue))+(2*(1-alpha)*U)
    
#Similar to AX=B here we have left_matrix*temp=coeff_matrix
    temp= np.linalg.inv(left_matrix).dot(coeff_matrix)
    for k in range(1,length_nodes+1):
            grid[i-1,k]=round(temp[k-1,0],2)                #after calculating temp it is updated in the grid
    
# print(grid)


# ----------------------------------------------------------------------------------------------------------------------------
#                                                          OUTPUT
# ----------------------------------------------------------------------------------------------------------------------------

print('Temperature of rod at the end of simulation time: ')
for i in range(n):
    print(grid[0,i],end=" ")

#Output in graph form
#X axis-length ; Y axis-Temperature
#Graph shows 5 plots with time difference=simulation time/5


if time_nodes>5:          
    X=[]
    Y=[]
    d=int(time_nodes/5)                                                                 #divides total time into 5

    time=t 

    for i in range(0,m,d):

        X.append(0)
        Y.append(t1)
        
        for j in range(1,length_nodes+1):
            
            Y.append(grid[i,j])
            X.append(j*dx)

        Y.append(t2)
        X.append((length_nodes+1)*dx) 

        s="at t ="+str(time)
        plt.plot(X,Y,label=s) 

        time=round((time-(dt*d)),2)
        X.clear()
        Y.clear()
    plt.legend()
    plt.show()