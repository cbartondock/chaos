from optparse import OptionParser
import time
import operator
start = time.time()

parser = OptionParser()
parser.add_option('-f','--file',dest='filename')
(options, args) = parser.parse_args()

with open(options.filename) as f:
    lines = f.readlines();
    lines_dict= [(float(line.split("zeros: ")[1].strip()), line) for line in lines]
    lines_dict.sort(key= lambda x: x[0])
    lines = [e[1] for e in lines_dict]
    #print lines_dict
    with open(options.filename.split(".txt")[0]+"_sorted.txt","w") as out:
        for line in lines:
            out.write(line)
end = time.time()
print end-start




