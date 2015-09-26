name = "outputs/text_quasi_conv_t20000_g2000_xs0.00_ys0.00_xb6.28_yb6.28.txt"
f = open(name,"r")

f2 = open(name.split(".txt")[0]+"_fordas.txt","w")
for line in f.readlines():
    data = [float(item.split(": ")[1]) for item in line.split(', ')[2:]]
    f2.write("{0} {1} {2}\n".format(*data))




