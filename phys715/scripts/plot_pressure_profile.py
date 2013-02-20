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

def plot_p(data, label):
    p = []
    r = []
    for point in data:
        r.append(point['r'])
        p.append(-1.0*point['p'])
    print "Mean Pressure: ", label, ":", np.mean(p)
    plt.plot(r, p, 'o', label=label)

if __name__ == "__main__":
    radius = 0.05
    empty_csv = csv.DictReader(open('../data/empty.csv', 'r'))
    cylinder_csv = csv.DictReader(open('../data/cylinder.csv', 'r'))
    NACA5440_csv = csv.DictReader(open('../data/NACA5440.csv', 'r'))
    NACA5420_csv = csv.DictReader(open('../data/NACA5420.csv', 'r'))
    empty_ann = get_points_annulus(empty_csv, radius)
    cylinder_ann = get_points_annulus(cylinder_csv, radius)
    NACA5440_ann = get_points_annulus(NACA5440_csv, radius)
    NACA5420_ann = get_points_annulus(NACA5420_csv, radius)
    plot_p(empty_ann, "Empty Box")
    plot_p(NACA5440_ann, "NACA5440")
    plot_p(NACA5420_ann, "NACA5420")
    plt.xlabel("Radius (m)")
    plt.ylabel("Flow Pressure (Pa)")
    plt.legend(loc='best')
    plt.show()
