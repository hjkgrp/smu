# generate exhaustive output list for mono-ligands
dictMo = {}
scoredDictMo = {}
wishedDictMo = {}

for elem in elemList.keys():
    for charge in chargeList:
        l = ligand(elem, 0, 0)
        l.setCharge(charge)
        for h in hList:
            if h:
                l.addH()
            l.getSmiles()
            dictMo[l.SMILES] = [l.numE,l.numValE,l.numLP,l.charge,l.numberOfHs,\
                                l.testOctetRule(),l.testValenceShell(),l.score()]
            
for name, line in dictMo.items():
    octetScore = line[7]
    charge = line[3]
    numberOfHs = line[4]
    numValE = line[1]
    closedShell = int(not((numValE)%2))
    # Charge: +1 >= charge >= -3
    # Sterics: 4 >= Number of H at Coordinating Atom (CA)
    # Closed Shell only
    if 1 >= charge and charge >= -3 and 4 >= numberOfHs and closedShell == 1: 
        
        # Charge score
        if charge == 1:
            scoreCharge = 0
        if charge <= 0 and charge >= -2:
            scoreCharge = 3
        elif charge == -3:
            scoreCharge = 0

        # CA Sterics Score
        if numberOfHs == 4:
            scoreCa = 0 
        else:
            scoreCa = 3
    
        # Total score
        score = scoreCharge + octetScore + scoreCa
        
        # Dict with only scored ligands
        dictMo[name] = line + [scoreCharge] + [octetScore] + [scoreCa] + [score] 
        scoredDictMo[name] = line + [scoreCharge] + [octetScore] + [scoreCa] + [score] 
        
    else:
        score = 0.0
        dictMo[name] = line + [0] + [0] + [0] + [score]

print("Name: Charge + octet + CA = Score")

histCharge = list() # lists of charge, octet difference, occupancy of CA and score to see histograms
histOctet = list()
histCa = list()
histScore = list()
histVE = list()
threshold = 7 # define threshold for wishedDictMo
for name, props in scoredDictMo.items():        
    histCharge.append(props[-4])
    histOctet.append(props[-3])
    histCa.append(props[-2])
    histScore.append(props[-1])
    histVE.append(props[1])
    
    # Populate the wishedDictMo
    if props[-1] > threshold: 
        wishedDictMo[name] = props
    
    # Evaluate compounds from Spectrochemical Series
    for i in range(0, len(scSeriesMo)):
        if name == scSeriesMo[i]:
            print(name + ': ' +str(props[-4])+" + "+str(props[-3])+" + "+str(props[-2])+' = '+str(props[-1])) 

print('======')
print("All monoatoms: " + str(len(dictMo)))  
print("All scored monoatoms: " + str(len(scoredDictMo)))
print("All wished for monoatoms (>" + str(threshold) + "): " + str(len(wishedDictMo)))
print('======')

plt.xlabel('Score')
plt.ylabel('Number of Wished for Mono-Ligands')
# plt.yscale('log', nonposy='clip')
plt.title('Histogram of Wished for Mono-Ligands ['+ components + "]")
plt.hist(histScore)
# plt.savefig('monoDistr' + components + ".pdf", bbox_inches='tight')
plt.show()


plt.hist2d(histVE, histScore, bins=4)
plt.xlabel('VE')
plt.ylabel('score')
plt.colorbar()
plt.show()
    
