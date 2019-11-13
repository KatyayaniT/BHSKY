import numpy as np

f1 = open("star.dat", "r+")
star_data = f1.readlines()
f1.close()

test = [[]]
for i in range(0, 9):
    test.insert(i, str(star_data[i+1]))

j = (np.float(test[5][20:29]))

while j <= 360.0:
    f1 = open("star.dat", "r+")

    star_data = f1.readlines()

    test = (str(star_data[1]))
    test1 = (str(star_data[2]))
    test2 = (str(star_data[3]))

    Rs = (np.float(test[20:29]))
    Rad = (np.float(test1[20:29]))
    Dobs = (np.float(test2[20:29]))

    f1.close()

    #exec(open("bhsky_table 3.py").read())
    exec(open("bhsky_image3.py").read())
    # exec(open("camera.py").read())

    star_data[6] = star_data[6].replace(str(j), str(j + 1.0))
    
    open('star.dat', "w").writelines(star_data)

    j = j + 1.0

