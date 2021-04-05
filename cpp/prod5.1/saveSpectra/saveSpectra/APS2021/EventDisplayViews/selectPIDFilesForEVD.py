import numpy as np

#CohArray = np.loadtxt('AllCCCoherentEvents.txt', dtype = str)
#CohArray = np.loadtxt('AllCCDISEvents.txt', dtype = str)
CohArray = np.loadtxt('AllCCRESEvents.txt', dtype = str)
#print(CohArray.shape)
#print("000"+str(CohArray[:,0])+"_s")
prongVertexE = []
keys = []

for row in CohArray:
    key = "000"+row[0]+"_s{:02d}".format(int(row[1]))+"_c{:03d}".format(int(row[2]))
    info = "Event: %s Slice: %s E_nu: %s VTX_z: %s, prong3D: %s"%(row[3], row[4], row[5],row[8],row[9])
    keyPack = []
    keyPack.append(key)
    keyPack.append(info)
    #print(info)
    #print(row[9],row[10])
    prongVertexE.append(float(row[9]))
    keys.append(keyPack)
    if (float(row[9]) > 0.0008 and float(row[9]) < 0.001):
        print(key)
        
#print(keys)

prongVertexE.sort()
fileList = []
with open('cachedFiles.txt','r') as file:
    lines = file.readlines()
    for line in lines:
        fileList.append(line.rstrip("\n"))
        #print(line[-53:70])

    file.close()
#print(fileList)

for key in keys:
    for file in fileList:
        #print(key[0],file)
        if (key[0] in file):
            print(file, key[1])

#print(prongVertexE)
    
    #print(key)

# with open ('AllCCCoherentEvents.txt', 'r') as file:
    

#     lines = file.readlines()
#     for line in lines:
#         print(line)
        