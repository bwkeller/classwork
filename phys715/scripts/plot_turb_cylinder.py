#!/usr/bin/python
import csv
import numpy as np
import matplotlib.pyplot as plt

def fmtpoint(point):
    p = {'x':float(point['Points:0']), 'y':float(point['Points:1']),
    'z':float(point['Points:2']), 'vx':float(point['Velocity:0']), 
    'vy':float(point['Velocity:1']), 'vz':float(point['Velocity:2']), 
    'te':float(point['TurbEner']), 'p':float(point['Pressure'])}
    p['r'] = np.sqrt(p['x']*p['x'] + p['y']*p['y'])
    return p

def get_points_annulus(data, radius):
    points = []
    for point in data:
        p = fmtpoint(point)
        if p['r'] < radius:
            points.append(p)
    return points

def plot_turb(data, label):
    te = []
    r = []
    for point in data:
        r.append(point['r'])
        te.append(point['te'])
    print "Mean Turbulent Energy Loss: ", label, ":", np.mean(te)
    plt.semilogy(r, te, 'o', label=label)

if __name__ == "__main__":
    radius = 0.05
    cylinder_csv = csv.DictReader(open('../data/cylinder.csv', 'r'))
    cylinder_ann = get_points_annulus(cylinder_csv, radius)
    plot_turb(cylinder_ann, "Cylinder")
    plt.xlabel("Radius (m)")
    plt.ylabel("Energy Lost to Turbulence (J)")
    plt.show()
