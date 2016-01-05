from optparse import OptionParser
import time
import operator
start = time.time()

parser = OptionParser()
parser.add_option('-f','--file',dest='filename')
(options, args) = parser.parse_args()

with open(options.filename) as f:
    lines = f.readlines();
    lines_dict= [(int(line.split(":")[1].split(',')[0].strip()),int(line.split(":")[2].split(",")[0].strip()),float(line.split(":")[5].split(",")[0].strip()), line) for line in lines]
    lines_dict=filter(lambda x: x[0]+x[1]<500 and x[0]+x[1]>100 and x[2]<8.,lines_dict)
    lines = [e[3] for e in lines_dict]
    #print lines_dict
    with open(options.filename.split(".txt")[0]+"_purged.txt","w") as out:
        for line in lines:
            out.write(line)
end = time.time()
print end-start




