#!/usr/bin/python3.5

def writeConfigFile(det1_thickness_um, det2_thickness_um, distance_between_det_mm,
		    window_thickness_um, window_gap_mm):
    with open('../src/detector_config.txt','w') as f:
        f.write(str(det1_thickness_um) + '\n')
        f.write(str(det2_thickness_um) + '\n')
        f.write(str(distance_between_det_mm) + '\n')
        f.write(str(window_thickness_um) + '\n')
        f.write(str(window_gap_mm) + '\n')
