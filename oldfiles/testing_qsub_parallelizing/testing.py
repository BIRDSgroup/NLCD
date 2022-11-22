import sys
#First argument for the number of runs, second for the initial run
j=int(sys.argv[2])
n=int(sys.argv[1])
f=open("tesr.txt","a")
for i in range(j,j+n): 
  print(sys.argv[1]," ",sys.argv[2])
  f.write(str(i)+","+str(n)+"\n")
#f.close()
