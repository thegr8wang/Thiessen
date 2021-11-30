import math

import matplotlib.pyplot as plt
import numpy as np
import xlrd as xlrd
from scipy.spatial import KDTree
from shapely.geometry import Polygon, Point, MultiPoint


def xls2points(filepath, sheet=0, start_col=0, end_col=None):
    """
    Read a xls by row and write an numpy.ndarray.
    :return ndarray
    """
    points = []
    data = xlrd.open_workbook(filepath)
    table = data.sheets()[sheet]
    for i in range(table.nrows):
        points.append(table.row_values(i, start_col, end_col))
    return np.array(points)


def points2buffer_points(points):
    """

    :param points:
    :return: ndarray
    """
    radius = 2.0  # radius单位为度，0.00001约等于1米
    convex_buffer = MultiPoint(points).convex_hull.buffer(radius)
    return np.concatenate((points, convex_buffer.exterior.coords), axis=0)


def txt2data(filepath, start_line=0, end_line=None, spilt=","):
    """ Read 'export_output.txt' and all txt like this. """
    data = []
    with open(filepath, 'r', encoding='utf-8') as f:
        f_lines = f.readlines()
        if end_line is None:
            end_line = len(f_lines)
        for line in f_lines[start_line:end_line]:
            line_list = line.strip('\n').split(spilt)
            point = [float(i) for i in line_list]
            data.append(point)
    return np.array(data)


def remove_duplicates(array):
    """

    :param array: list
    :return: ndarray
    """
    result = [array[0]]
    for i in array:
        if i not in result:
            result.append(i)
    return np.array(result)


def boundary_points(sites, boundary):
    """

    :param sites:
    :param boundary:
    :return: list
    """
    kt = KDTree(sites)
    result = [[] for _ in range(len(sites))]
    for i in boundary:
        (d, p) = kt.query(i)
        result[p].append(i)
    return result


def sort_with_insert(arr):
    result = []
    for a in arr:
        if not result:
            result.append(arr[0])
        else:
            for i in range(len(result)):
                if a > result[i]:
                    result.insert(i, a)
                    break
            if a not in result:
                result.append(a)
    return result


def sort_with_slope(arr, site):
    result = []
    slopes = []
    inf_yu = []
    inf_yd = []

    def sort_reserve_order(array1, array2):
        if not array1:
            array1.append(array2)
        else:
            for index in range(len(array1)):
                if array2[1] > array1[index][1]:
                    array1.insert(index, array2)
                    break
            # if array2 not in array1:
            #     array1.append(array2)
            flag = False
            for ar in array1:
                if (np.array(ar) == np.array(array2)).all():
                    flag = True
                    break
            if not flag:
                array1.append(array2)
        return array1

    for a in arr:
        if a[0] == site[0]:
            if a[1] > site[1]:
                inf_yu = sort_reserve_order(inf_yu, a)
            else:
                inf_yd = sort_reserve_order(inf_yd, a)
        else:
            k = (a[1] - site[1]) / (a[0] - site[0])
            if not slopes:
                slopes.append(k)
                result.append(a)
            else:
                for i in range(len(slopes)):
                    if k > slopes[i]:
                        slopes.insert(i, k)
                        result.insert(i, a)
                        break
                if k not in slopes:
                    slopes.append(k)
                    result.append(a)
    for yu in inf_yu:
        result.append(yu)
    for yd in range(len(inf_yd) - 1, -1, -1):
        result.insert(0, inf_yd[yd])
    return result


def screen_sort_with_slope(sites, vertices, points):
    epsilon = 10 ** (-1)
    polygon = Polygon(polygon_box_bound(points))
    polygon_vertices = [[] for _ in sites]
    sort_ver = [[] for _ in sites]
    for i in range(len(sites)):
        if points[i]:
            site = sites[i]
            slopes = []
            max_point_x = points[i][0]
            for point in points[i]:
                if point[0] > max_point_x[0]:
                    max_point_x = point
                if point[0] == site[0]:
                    k = math.exp(site[0])
                else:
                    k = (point[1] - site[1]) / (point[0] - site[0])
                if k not in slopes:
                    slopes.append(k)
                    polygon_vertices[i].append(point)
            max_vertice_x = max_point_x
            for vertice in vertices[i]:
                if polygon.intersects(Point(vertice)):
                    polygon_vertices[i].append(vertice)
                    if vertice[0] > max_vertice_x[0]:
                        max_vertice_x = vertice
                    if vertice[0] == site[0]:
                        kv = math.exp(site[0])
                    else:
                        kv = (vertice[1] - site[1]) / (vertice[0] - site[0])
                    flag = False
                    for slope in slopes:
                        if kv + epsilon >= slope >= kv - epsilon:
                            flag = True
                            break
                    if not flag:
                        slopes.append(kv)
                        polygon_vertices[i].append(vertice)
            max_x = max_point_x if max_point_x[0] > max_vertice_x[0] else max_vertice_x
            for pv in range(len(polygon_vertices[i])):
                if (np.array(polygon_vertices[i][pv]) == np.array(max_x)).all():
                    polygon_vertices[i].pop(pv)
                    break
            sort_ver[i] = sort_with_slope(polygon_vertices[i], max_x)
            sort_ver[i].append(max_x)
            sort_ver[i].insert(0, max_x)
            # print(slopes)
            # print(max_x)
        else:
            sort_ver[i] = vertices[i]
    return sort_ver


def screen_sort(sites, vertices, points):
    polygon_vertices = [[] for _ in sites]
    sort_ver = [[] for _ in sites]
    polygon = Polygon(polygon_box_bound(points))
    for i in range(len(sites)):
        if points[i]:
            polygon_vertices[i] = points[i]
            max_point_x = points[i][0]
            for point in points[i]:
                if point[0] > max_point_x[0]:
                    max_point_x = point
            max_vertice_x = max_point_x
            for vertice in vertices[i]:
                if polygon.intersects(Point(vertice)):
                    polygon_vertices[i].append(vertice)
                    if vertice[0] > max_vertice_x[0]:
                        max_vertice_x = vertice
            max_x = max_point_x if max_point_x[0] > max_vertice_x[0] else max_vertice_x
            for pv in range(len(polygon_vertices[i])):
                if (np.array(polygon_vertices[i][pv]) == np.array(max_x)).all():
                    polygon_vertices[i].pop(pv)
                    break
            sort_ver[i] = sort_with_slope(polygon_vertices[i], max_x)
            sort_ver[i].append(max_x)
            sort_ver[i].insert(0, max_x)
        else:
            sort_ver[i] = vertices[i]
    return sort_ver


def polygon_box_bound(points):
    [max_x, max_y, min_x, min_y] = [None for _ in range(4)]

    for point in points:
        if point:
            if max_x is None:
                max_x = point[0][0]
            max_x = max(max(p[0] for p in point), max_x)
            if min_x is None:
                min_x = point[0][0]
            min_x = min(min(p[0] for p in point), min_x)
            if max_y is None:
                max_y = point[0][1]
            max_y = max(max(p[1] for p in point), max_y)
            if min_y is None:
                min_y = point[0][1]
            min_y = min(min(p[1] for p in point), min_y)
    xy = [max_x, max_y, min_x, min_y]
    box_bound = [[xy[0], xy[1]], [xy[0], xy[3]], [xy[2], xy[3]], [xy[2], xy[1]], [xy[0], xy[1]]]
    return box_bound


def plot_lines(arr):
    for a in arr:
        plt.plot([v[0] for v in a], [v[1] for v in a], 'o', linewidth=1)
    plt.show()


def cal_area(arr):
    areas = []
    for a in arr:
        polygon = Polygon(a)
        areas.append(polygon.area)
    return areas
