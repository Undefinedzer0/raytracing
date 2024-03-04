import numpy as np
import scene
from scipy.spatial.transform import Rotation as R
import time
from multiprocessing import Pool


def IntersectRaySphere(O, D, sphere):
    C = sphere.center
    r = sphere.radius
    oc = O - np.array(C)
    k1 = np.dot(D, D)
    k2 = 2 * np.dot(oc, D)
    k3 = np.dot(oc, oc) - r * r
    dsc = k2 * k2 - 4 * k1 * k3
    if dsc < 0:
        return inf, inf
    t1 = (-k2 + dsc ** 0.5) / (2 * k1)
    t2 = (-k2 - dsc ** 0.5) / (2 * k1)
    return t1, t2


def length(m):
    return sum(i ** 2 for i in m) ** 0.5


def ComputeLighting(P, N, V, s):
    i = 0.0
    for light in lights:
        if light.typ == "ambient":
            i += light.intensity
        else:
            if light.typ == "point":
                L = np.array(light.position) - P
                t_max = 1
            else:
                L = np.array(light.direction)
                t_max = inf
            n_dot_l = np.dot(N, L)
            shadow_sphere, shadow_t = ClosestIntersection(P, L, 0.001, t_max)
            if shadow_sphere != None:
                continue
            if n_dot_l > 0:
                i += light.intensity * n_dot_l / (length(N) * length(L))
            if s != -1:
                R = 2 * N * np.dot(N, L) - L
                r_dot_v = np.dot(R, V)
                if r_dot_v > 0:
                    i += light.intensity * pow(r_dot_v / (length(R) * length(V)), s)
    return i


def ClosestIntersection(O, D, t_min, t_max):
    closest_t = inf
    closest_sphere = None
    for sphere in spheres:
        t1, t2 = IntersectRaySphere(O, D, sphere)
        if t_min <= t1 <= t_max and t1 < closest_t:
            closest_t = t1
            closest_sphere = sphere
        if t_min <= t2 <= t_max and t2 < closest_t:
            closest_t = t2
            closest_sphere = sphere
    return closest_sphere, closest_t


def TraceRay(O, D, t_min, t_max, depth):
    closest_sphere, closest_t = ClosestIntersection(O, D, t_min, t_max)
    if closest_sphere == None:
        return background_color
    P = O + closest_t * D  # вычисление пересечения
    N = P - np.array(
        closest_sphere.center
    )  # вычисление нормали сферы в точке пересечения
    N = N / length(N)
    local_color = np.array(closest_sphere.color) * ComputeLighting(
        P, N, -D, closest_sphere.specular
    )
    r = closest_sphere.reflectance
    if depth <= 0 or r <= 0:
        return local_color
    R = ReflectRay(-D, N)
    reflected_color = np.array(TraceRay(P, R, 0.001, inf, depth - 1))
    return (1 - r) * local_color + reflected_color * r


def ReflectRay(R, N):
    return 2 * N * np.dot(N, R) - R


def draw_pixel(xy):
    x, y = xy
    D = np.array([x * Vw / Cw, y * Vh / Ch, d])
    D = R.from_euler("xyz", rotation, degrees=True).apply(D)
    color = TraceRay(O, D, 1, inf, recur_depth)
    return [[x + Cw // 2, Ch // 2 - y - 1], [max(min(255, i), 0) for i in color]]


def draw_screen():
    for x in range(-Cw // 2, Cw // 2, pool_size):
        print(x)
        for y in range(Ch // 2 - 1, -Ch // 2 - 1, -pool_size):
            pygame.event.get()
            for i in p.map(draw_pixel, [[x + xx, y - yy] for xx in range(pool_size) for yy in range(pool_size)]):
                screen.set_at(*i)
            pygame.display.flip()

def draw_screen_but_cooler():
    n = set()
    for y in range(Ch // 2 - 1, -Ch // 2 - 1, -1):
        for x in range(-Cw // 2, Cw // 2):
            n.add((x, y))
    n = list(n)
    for j in range(0, Ch * Cw, pool_size**2):
        pygame.event.get()
        for i in p.map(draw_pixel, [n[j + h] for h in range(pool_size**2)]):
            screen.set_at(*i)
        pygame.display.flip()


res = []
spheres = scene.spheres()
lights = scene.lights()
run = True
Cw = 500
Ch = 500
Vw = 1
Vh = 1
d = 1
inf = 2 ** 20
recur_depth = 3
background_color = [255, 255, 255]
rotation = [0, 0, 0]
O = [0, 0, 0]
pool_size=20
if __name__ == "__main__":
    import pygame
    
    p = Pool(pool_size)
    pygame.init()
    clock = pygame.time.Clock()
    screen = pygame.display.set_mode((Cw, Ch))
    pygame.display.set_caption("fast rt")
    start = time.time()
    cool = False
    if cool:
        draw_screen_but_cooler()
    else:
        draw_screen()
    print(time.time() - start)
    pygame.display.flip()
    while run:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
        clock.tick(60)
    pygame.quit()
