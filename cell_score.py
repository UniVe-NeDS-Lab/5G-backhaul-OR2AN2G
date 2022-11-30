#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 11:13:56 2022

@author: leonardo
"""

from numpy.random import uniform
from shapely.geometry import Point, Polygon, MultiPoint
import matplotlib.pyplot as plt
from collections import Counter, defaultdict



def set_scenario(radius, coords={}, size=0, n=0):
    """ Creates a scenario with randomly placed nodes and generates
    shapely circles to represent coverage. Can receive the list of
    positions as a parameter """
    if size and size < 2*radius:
        print('ERROR: area side must be at least twice the coverage radius')
        exit()
    circles = []
    if not coords:
        for i in range(1,n+1):
            x = uniform(0+radius,size-radius)
            y = uniform(0+radius,size-radius)
            coords[i] = (x,y)
    for k,v in coords.items():
        x,y = v
        p = Point(x, y).buffer(radius)
        circles.append((k,x,y,radius,p))
    return circles

def make_square(x,y,side=10):
    """ Just make a square from the center point. Useful to plot """
    return Polygon([(x-side/2, y-side/2), (x-side/2, y+side/2), 
                   (x+side/2, y+side/2), (x+side/2, y-side/2)])

def points_in_circle(x, y, radius, step=10, plot_squares=False):
    """ extract all coorsdinates of points inside a circle, sampling 
    every step meters. Returns points, convex hull and squares if needed. """
    samples = []
    squares = []
    min_x = int(((x-radius)//step)*step)-step
    max_x = int(((x+radius)//step)*step)+step

    for xi in range(min_x, max_x, step):
        delta = (radius**2 - (xi-x)**2)
        if delta < 0:
            continue
        max_y_f = y + delta**0.5
        min_y_f = y - delta**0.5
        min_y = int((min_y_f//step)*step)
        max_y = int((max_y_f//step)*step)+step
        for yi in range(min_y, max_y, step):
            if plot_squares:
                squares.append(make_square(xi,yi,step))
            samples.append((xi,yi))
            
    return samples, MultiPoint(samples).convex_hull, squares
            
            
def sample(circles,step=10):
    """ Sample all the area, and return a weight for every circle, and a map
    point -> circle """ 
    points_multiplicity = Counter()
    circle_to_points = defaultdict(list)
    for i,x,y,r,_ in circles:
        samples, union, _= points_in_circle(x,y,r,step)
        for point in samples:
            points_multiplicity[point] += 1
            circle_to_points[i].append(point)
    circle_score = {}
    max_score = 0
    for c in circles:
        idx = c[0]
        score = sum([1/points_multiplicity[p] for p in circle_to_points[idx]])
        circle_score[idx] = int(score)
        max_score = max(max_score, circle_score[idx])
    # the weight of each gNB must be an integer, and we normalize to 
    # 1000. Note that we need a capacity per node higher than that number.
    for k,v in circle_score.items():
        circle_score[k] = int(1000*v/max_score)
    return circle_score, points_multiplicity
            

def plot(circles,step, mult={}):
    """ just draw the area, for debugging purpuses """
    for i,x,y,r,p in circles:
        x_shape,y_shape = p.exterior.xy
        plt.plot(x_shape,y_shape)
        samples, union, squares = points_in_circle(x,y,r,step)
        for s in squares:
            plt.plot(*s.exterior.xy)
        #plt.plot(*union.exterior.xy)
    if mult:
        x_coord = []
        y_coord = []
        color = []
        m = max(mult.values())
        for x,y in mult:
            x_coord.append(x)
            y_coord.append(y)
            color.append(mult[(x,y)]/m)
        plt.scatter(x_coord,y_coord, alpha=0.5, c=color, cmap='Reds')
    plt.show()
    
    
if __name__ == "__main__":
    from manage_graphs import make_poisson_graph
    g, coords, seed = make_poisson_graph(40)
    radius = 100
    side=10
    print(coords, radius)
    c = set_scenario(coords=coords, radius=radius)
    #c = set_scenario(radius=100,size=1000,n=1)
    _, mult = sample(c,side)
    plot(c,side,mult)
