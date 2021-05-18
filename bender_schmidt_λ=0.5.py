import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# INPUT
# ----------------------------------------------------------

k=float(input('Enter diffusion coefficient(in cm^2/s): '))
x=float(input('Enter length of rod(in cm): '))
t=float(input('Enter simulation time(in s): '))

#node sizes
dx=1
dt=0.5*k

#calculating no. of nodes based on node size
length_nodes=int(x/dx)-1
time_nodes=int(t/dt)
alpha= 0.5

#calculating number of rows and columns for the grid which will show temperature at a node n at time t

grid_rows=time_nodes+1
grid_columns=length_nodes+2
grid=np.zeros((grid_rows,grid_columns))

#initialising the grid with initial and boundary conditions

print('Enter boundary conditions(in celcius')
T1=int(input('Enter surface 1 temperature: '))
T2=int(input('Enter surface 2 temperature: '))

for i in range(grid_rows):
    grid[i,0]=T1
    grid[i,length_nodes+1]=T2
t_initial=int(input('Enter initial rod temperature: '))
for i in range(1,length_nodes+1):
    grid[time_nodes,i]=t_initial


# ----------------------------------------------------------
#CALCULATION
# ----------------------------------------------------------

#first loop increments row of the grid(i.e)Time
#second loop increments column of the grid(i.e)length

for i in range(time_nodes,0,-1):
    for j in range(1,length_nodes+1):

        Uiw=grid[i,j-1]                             #temperature at node to the west
        Uie=grid[i,j+1]                             #temperature at node to the east
       
        Ujn=alpha*(Uie+Uiw)                          #calculates temperature at current node at time t+1 using finite difference formula
        grid[i-1,j]=round(Ujn,4)                     #calculated value is updated to the grid

# print(grid)



# ----------------------------------------------------------
#OUTPUT
# ----------------------------------------------------------

print('Temperature of rod at the end of simulation time: ')
for i in range(grid_columns):
    print(grid[0,i],end=" ")

#Output in graph form
#X axis-length ; Y axis-Temperature
#Graph shows 5 plots with time difference=simulation time/5

if time_nodes>5:
    X=[]
    Y=[]
    d=int(time_nodes/5)                                                 #divides total time into 5
    time=t 
    for i in range(0,grid_rows,d):
        
        X.append(0)
        Y.append(T1)
        
        for j in range(1,length_nodes+1):
    
            Y.append(grid[i,j])
            X.append(j*dx)

        Y.append(T2)
        X.append((length_nodes+1)*dx) 

        s="at t ="+str(time)
        plt.plot(X,Y,label=s) 
        time=round((time-(dt*d)),2)
        X.clear()
        Y.clear()
    plt.legend()
    plt.show()







