"""
Calculate the weight of the watershed surface
controlled by the rainfall station
incorporating Thiessen Diagram
"""
import numpy as np
import timeit as tt
from scipy.spatial import Voronoi

import thiessen as tn

""" Clazz Form """


class ThiessenDiagram:
    def __init__(self, sites, boundary):
        self.area_power = None
        self.sites = sites
        self.boundary = boundary

    def thiessen(self, screen_with_slopes=True):
        sites = self.sites
        boundary = self.boundary
        buffer_sites = tn.points2buffer_points(sites)
        vor = Voronoi(buffer_sites, furthest_site=False, incremental=True, qhull_options=None)
        vertices = vor.vertices
        point_region = vor.point_region
        regions = vor.regions
        polygon_vertices = [[] for _ in sites]  # Type:list
        for index, region in enumerate(regions):
            if len(region) == 0:
                continue
            if -1 in region:
                continue
            pv = []
            for r in region:
                pv.append(vertices[r])
            pv.append(vertices[region[0]])
            polygon_vertices[np.where(vor.point_region == index)[0][0]] = pv
        boundary = tn.remove_duplicates(boundary.tolist())
        polygon_points = tn.boundary_points(sites, boundary)
        if screen_with_slopes:
            result = tn.screen_sort_with_slope(sites, polygon_vertices, polygon_points)
        else:
            result = tn.screen_sort(sites, polygon_vertices, polygon_points)
        area = tn.cal_area(result)
        self.area_power = [i / sum(area) for i in area]


""" Function Form """


def thiseen_diagram(sites, boundary, screen_with_slopes=True):
    buffer_sites = tn.points2buffer_points(sites)
    vor = Voronoi(buffer_sites, furthest_site=False, incremental=True, qhull_options=None)
    vertices = vor.vertices
    regions = vor.regions
    polygon_vertices = [[] for _ in sites]  # Type:list
    for index, region in enumerate(regions):
        if len(region) == 0:
            continue
        if -1 in region:
            continue
        pv = []
        for r in region:
            pv.append(vertices[r])
        pv.append(vertices[region[0]])
        polygon_vertices[np.where(vor.point_region == index)[0][0]] = pv
    boundary = tn.remove_duplicates(boundary.tolist())
    polygon_points = tn.boundary_points(sites, boundary)
    if screen_with_slopes:
        result = tn.screen_sort_with_slope(sites, polygon_vertices, polygon_points)
    else:
        result = tn.screen_sort(sites, polygon_vertices, polygon_points)
    area = tn.cal_area(result)
    return [i / sum(area) for i in area]


if __name__ == '__main__':
    start_time = tt.default_timer()
    path_sites = r'sites.txt'
    sites = tn.txt2data(path_sites, start_line=1)
    path_boundary_points = r'boundary_points.txt'
    boundary = tn.txt2data(path_boundary_points, start_line=1)
    ts = ThiessenDiagram(sites, boundary)
    ts.thiessen()
    ap = ts.area_power
    end_time = tt.default_timer()
    print(ap)
    print("Time:" + '\t' + str(end_time - start_time))
    print(thiseen_diagram(sites, boundary))
    