import GeneParameter as GP
import gene as gene
import os

path = './data/'
dataList = []
input = []
output = []
for dirPath, dirNames, fileNames in os.walk(path):
    for file in fileNames:
        file = os.path.join(dirPath, file)
        f = open(file, 'r')
        print(file)
        for line in f.readlines():
            line_data = line.split()
            mylist = []
            if type(eval(line_data[0])) == float:
                mylist.append(eval(line_data[0]))
                mylist.append(eval(line_data[1]))
                mylist.append(eval(line_data[2]))
                input.append(mylist)
                output.append(eval(line_data[3]))
            else:
                print(line_data)

for i in range(0, len(output)):    # Normalization(0~1)
    output[i] = (output[i] + 40) / 80

fError_ori = 1e9
fError_now = 1e9

while fError_now > 0.1:
    if os.path.isfile("./bestGP.txt"):
        readFile = open("./bestGP.txt", 'r')
        gene_ori = readFile.read()  # 16 parameters
        readFile.close()
        g = gene.Gene()
        g.setGene(gene_ori)
        g.calculateFitness(input, output)
        fError_ori = g.f
        gp = GP.GeneParameter(bestGene = g)
    else:
        gp = GP.GeneParameter()

    g = gp.geneIteration(input, output)
    g.calculateFitness(input, output)

    fError_now = g.f

    if fError_now < fError_ori:
        writeFile = open("./bestGP.txt", 'w')
        for i in range(len(g.DNA)):
            writeFile.write(str(g.DNA[i]))
            writeFile.write('     ')
        writeFile.close()
        print("file has be written")
    print("--------------------------------------")