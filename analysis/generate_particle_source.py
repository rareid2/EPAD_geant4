import random as rand

rand.seed()	# starts the Mersenne Twister rand engine

# Constant offset upstream of detector
y_offset = -5	# cm
num_sources = 100

# Maximum X and Z coordinate ranges that can be randomly chosen
x_range = 5
z_range = 5

counter = 1

with open('run_rand_sources.mac', 'w') as f:
    # init lines to set up simulation
    f.write('/run/initialize \n') 
    f.write('/control/verbose 0 \n')
    f.write('/run/verbose 0 \n')
    f.write('/event/verbose 0 \n')
    f.write('/tracking/verbose 0 \n')
    
    for i in range(0, num_sources):
    
        x_pos = round(rand.uniform(-x_range/2, x_range/2), 3)
        z_pos = round(rand.uniform(-z_range/2, z_range/2), 3)
        position_string = str(x_pos) + ' ' + str(y_offset) + ' ' + str(z_pos) + ' cm'
        
        x_dir = round(rand.uniform(-1, 1), 3)
        y_dir = 1 # maintains the general direction towards the detector
        z_dir = round(rand.uniform(-1, 1), 3)
        direction_string = str(x_dir) + ' ' + str(y_dir) + ' ' + str(z_dir)
	
        f.write('\n# Particle source ' + str(counter) + '\n')
        f.write('/gps/particle e-\n')
        f.write('/gps/ene/min 100. keV \n')
        f.write('/gps/ene/max 10. MeV \n')
        f.write('/gps/position ' + position_string + '\n')
        f.write('/gps/direction ' + direction_string + '\n')
        f.write('/gps/pos/type Point \n')
        
        f.write('\n/run/beamOn 100 \n')

        counter = counter + 1

    
    f.close()

