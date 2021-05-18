import numpy as np
import matplotlib.pyplot as plt

# INPUT

k=0.835
x=10
nx=4
t=50
nt=500


alpha= 0.0208
diag=2*(1+alpha)
below_diag=-alpha
above_diag=-alpha

m=nt+1
n=nx+2
grid=np.zeros((m,n))


t1=100
t2=50

for i in range(m):
    grid[i,0]=round(t1,2)
    grid[i,nx+1]=round(t2,0)
t_initial=0
for i in range(1,nx+1):
    grid[nt,i]=round(t_initial,2)

#print('initial conditions ',grid)
left_matrix=np.zeros((nx,nx))
for i in range(nx-1):
            left_matrix[i][i]=diag
            left_matrix[i][i+1]=above_diag
            left_matrix[i+1][i]=below_diag
left_matrix[nx-1][nx-1]=diag

print(left_matrix)

#PROCESSING
temp=np.zeros((nx,1))
d=np.zeros((nx,1))
for i in range(nt,1,-1):
    for j in range(1,nx+1,1):
        if j==1:
            U=grid[i,j]
            Uw=grid[i,j-1]
            Ue=grid[i,j+1]
            Unw=grid[i-1,j-1]
            Une=grid[i-1,j+1]
            d[j-1,0]=(alpha*(Uw+Unw+Ue))+(2*(1-alpha)*U)
            # print('1')

        elif j==(nx):
            U=grid[i,j]
            Uw=grid[i,j-1]
            Ue=grid[i,j+1]
            Une=grid[i-1,j+1]
            d[j-1,0]=(alpha*(Uw+Ue+Une))+(2*(1-alpha)*U)
            # print('2')
         

        else:
            U=grid[i,j]
            Uw=grid[i,j-1]
            Ue=grid[i,j+1]
            d[j-1,0]=(alpha*(Uw+Ue))+(2*(1-alpha)*U)
            # print('4')
    #print(d)
    temp= np.linalg.inv(left_matrix).dot(d)
    for k in range(1,nx+1):
            grid[i-1,k]=round(temp[k-1,0],2)
    # print(temp)
print(temp)







# for i in range(nt,0,-1):
#     for j in range(1,nx+1):
#         Ui=grid[i,j]
#         Uiw=grid[i,j-1]
#         Uie=grid[i,j+1]
#         #print(Uie,Ui,Uiw,'i',i,'j',j)
#         Ujn=Ui+(alpha*(Uie-2*Ui+Uiw))
#         grid[i-1,j]=round(Ujn,4)
#         print(Ujn)
#     print('next')

# d=int(t/5)
# for i in np.arange(0,t,d):
#     print('0')
#     X=[0]
#     Y=[100]
#     for j in range(1,n):
#         #print(i,j)
#         Y.append(grid[i,j])
#         X.append(j)
#     plt.plot(X,Y)
#     X.clear()
#     Y.clear()
# plt.show()


# d=int(nt/5)
# for i in np.arange(0,nt,d):
#     # print('0')
#     X=[0]
#     Y=[100]
#     for j in range(1,n):
#         #print(i,j)
#         Y.append(grid[i,j])
#         X.append(j)
#     print('X',X)
#     print('Y',Y)
#     plt.plot(X,Y)
#     X.clear()
#     Y.clear()
#     print(X)
# plt.show()


