import math

nTheta = 20
nL = 20

L = 600.
R = 300.

dTheta = math.pi / nTheta
dL = L / nL
for i in range(-nTheta // 2, nTheta // 2 + 1):
    for j in range(-nL // 2, nL // 2 + 1):
        theta = dTheta * i
        print("v",
              R * math.sin(theta),
              dL * j,
              R * math.cos(theta) - R * math.cos(40 / 180 * math.pi))


def index(i, j):
    return i * (nL + 1) + j


for i in range(nTheta):
    for j in range(nL):
        print("f", index(i, j) + 1, index(i, j+1) + 1, index(i+1, j+1) + 1, index(i+1, j) + 1)

fix_xz = [index(i, 0) for i in range(nTheta+1)] + [index(i, nL) for i in range(nTheta+1)]
fix_z = [index(0, i) for i in range(nL+1)] + [index(nTheta, i) for i in range(nL+1)]

fix_xz = list(map(str, fix_xz))
fix_z = list(map(str, fix_z))

print("# fixxz", ", ".join(fix_xz))
print("# fixz", ", ".join(fix_z))
print("# measure", index(nTheta // 2, nL // 2))
