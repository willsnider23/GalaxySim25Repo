import numpy
import visual
import linecache
import time
import math

#Variables
outputHead = 1  #Line number of the start of the next output section
N = 0           #Number of bodies
point = [] 
datarray = []   #2D array to store all position values, used to update body positions after each output section
time = 0        #Time of current output section
Nbod = []       #List of bodies used by visualpython to draw animation frames
outputNum = 0   #Counter tracking which output section is being read
r_half = 0      #Half-mass radius at start
escapee_ids = []

size = 750
origin = -size*3/4

#Visual#
scene = visual.display(title='Dsph',x=0,y=0,width=1200,height=1200,range=size,autoscale=1,userspin = True) #change x location to width/n
rod = visual.cylinder(pos=(origin,origin,0), axis=(0,1000,0), radius=10, color=visual.color.orange)
rod = visual.cylinder(pos=(origin,origin,0), axis=(1000,0,0), radius=10, color=visual.color.cyan)

L=visual.label(text = str(0) , pos=(0, -7000, 0), height = 20, color=visual.color.yellow)

#Identify Escapers
#with open("escapers.txt") as esc:
#    for line in esc:
#        data = line.split()
#        escapee_ids.append(float(data[0]))

#Main Program#
linenum = 1
line = " "
f = open("Run_2results.txt", "r")
# Until the end of the file, read line by line and process based on placement in output blocks
while line != "end":
    line = f.readline()
    visual.rate(10000)
    # print(line);
    values=line.split()
    # At the start of each code block, update the display time and the number of bodies N
    # then increment the outputHead tracker
    if linenum == outputHead:
        outputNum = outputNum + 1
        time = values[0]
        #print(time)
        N = int(values[1])
        if (outputNum == 1):
            r_half = float(values[2])
        L.text=time
        outputHead = outputHead + N + 1
    # Within each code block, for each star read the position coordinates (x, y, z)
    elif (float(values[0]) not in escapee_ids):
        r = float(values[1]) #math.sqrt(float(values[1])**2 + float(values[2])**2 + float(values[3])**2)
        v = float(values[4]) #math.sqrt(float(values[4])**2 + float(values[5])**2 + float(values[6])**2)
        point.append((r/2) + origin)      #(1000*math.log10(r/r_half) + origin)
        point.append(10*v + origin)
        # if linenum % 2 == 0:
        datarray.append(point)
        if (outputNum == 1):
            Nbod.append(visual.sphere(pos=(point[0],point[1],0),radius=20,color=(1,1,1),make_trail=0)) #, retain=300))
        point = []
    #Before the start of a new output section, take all bodies and update their positions based on datarray values
    if linenum == outputHead - 1:
        for i in range(len(Nbod)):
           #if (i in escapee_ids):
                Nbod[i].pos=(datarray[i][0],datarray[i][1],0)
                Nbod[i].make_trail=1
        datarray = []
    linenum = linenum + 1
