import numpy
import visual
import linecache
import time

#Variables
outputHead = 1  #Line number of the start of the next output section
N = 0           #Number of bodies
positions = []  #To store the position vector of each body [x, y, z]
datarray = []   #2D array to store all position values, used to update body positions after each output section
time = 0        #Time of current output section
prev_time = 0   #Prev output time
Nbod = []       #List of bodies used by visualpython to draw animation frames
outputNum = 0   #Counter tracking which output section is being read
bound = []
unbound = []
boundedness = []

#Visual#
scene = visual.display(title='Newton',x=0,y=0,width=800,height=1050,range=2500,autoscale=1,userspin = True,ambient=visual.color.gray(.4)) #change x location to width/n
#rod = visual.cylinder(pos=(0,0,-300), axis=(0,0,600), radius=1, color=visual.color.yellow)
rod = visual.cylinder(pos=(0,0,0), axis=(1500,0,0), radius=10, color=visual.color.yellow)
rod = visual.cylinder(pos=(0,0,0), axis=(0,1500,0), radius=10, color=visual.color.orange)
rod = visual.cylinder(pos=(0,0,0), axis=(0,0,1500), radius=10, color=visual.color.cyan)
ring = visual.ring(pos=(0,0,0), axis=(0,0,1), radius=20, thickness=2, color=visual.color.red)#
ring = visual.ring(pos=(0,0,0), axis=(0,0,1), radius=1000, thickness=4, color=visual.color.blue)
ring = visual.ring(pos=(0,0,0), axis=(0,0,1), radius=2000, thickness=4, color=visual.color.red)
ring = visual.ring(pos=(0,0,0), axis=(0,1,0), radius=1000, thickness=4, color=visual.color.blue)
ring = visual.ring(pos=(0,0,0), axis=(0,1,0), radius=2000, thickness=4, color=visual.color.red)
#ring = visual.ring(pos=(2000,0,0), axis=(0,0,1), radius=500, thickness=25, color=visual.color.white)
#ring = visual.ring(pos=(2000,0,0), axis=(0,1,0), radius=500, thickness=25, color=visual.color.white)
L=visual.label(text = ' ' , pos=(0, -3000, 0), height = 20, color=visual.color.yellow)

#Main Program#
linenum = 1
line = " "
f = open("results.txt", "r")    #Run_1
# Until the end of the file, read line by line and process based on placement in output blocks
while line != "end":
    line = f.readline()
    visual.rate(50000)
    values=line.split()
    # At the start of each output block, update the display time and the number of bodies, N
    # then increment the outputHead tracker
    if linenum == outputHead:
        outputNum = outputNum + 1
        if outputNum == 1:
            prev_time = time
        timetxt = values[0]
        time = float(timetxt)
        N = int(values[1])
        boundedness = ['b'] * N
        # read COM position and enter as first object in array
        #COM_vx = float(values[3])
        #COM_vy = float(values[4])
        #COM_vz = float(values[5])
        L.text=timetxt
        outputHead = outputHead + N + 1
    # Within each output block, for each star read the position coordinates (x, y, z)
    else:
        positions.append(float(values[1]))
        positions.append(float(values[2]))
        positions.append(float(values[3])) 
        datarray.append(positions)
        #if (values[7] == 'b' and float(values[0]) not in bound):
        #    if (float(values[0]) in unbound):
        #        unbound.remove(float(values[0]))
        #    bound.append(float(values[0]))
        #elif (values[7] == 'u' and float(values[0]) not in unbound):
        #    if (float(values[0]) in bound):
        #        unbound.remove(float(values[0]))
        #    unbound.append(float(values[0]))
        if outputNum == 1:
            boundedness[int(values[0])] = values[7]
            Nbod.append(visual.sphere(pos=(positions[0],positions[1],positions[2]),radius=10,color=(1,1,1),make_trail=1,retain=200))
        else:
            if (boundedness[int(values[0])] != values[7]):
                boundedness[int(values[0])] = values[7]
                if (values[7] == 'b'):
                    Nbod[int(values[0])].color=(1,1,1)
                elif (values[7] == 'u'):
                    Nbod[int(values[0])].color=(1,0,0)
        positions = []
    #Before the start of a new output section, take all bodies and update their positions based on datarray values
    if linenum == outputHead - 1:
        for i in range(N):
            Nbod[i].pos=(datarray[i][0],datarray[i][1],datarray[i][2])
        datarray = []
    linenum = linenum + 1
