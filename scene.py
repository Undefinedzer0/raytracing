class sphere:
    def __init__(self, cent, rad, col, spec, refl):
        self.center = cent
        self.radius = rad
        self.color = col
        self.specular = spec
        self.reflectance = refl
class light:
    def __init__(self, tp, ints, pos=None, dir=None):
        self.typ = tp
        self.intensity = ints
        if pos:
            self.position = pos
        if dir:
            self.direction = dir
r = []
with open('scene.txt') as f:
    r = f.read().split('\n')
def spheres():
    res = []
    for i in r:
        if 'sphere' in i:
             res.append(eval(i))
    return res
def lights():
    res = []
    for i in r:
        if 'light' in i:
             res.append(eval(i))
    return res