import math

nTheta = 240
nL = 400

L = 50.
R = 25.

dTheta = 80 / 360 * 2 * math.pi / nTheta
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

fixed = [index(i, 0) for i in range(nTheta+1)] + [index(i, nL) for i in range(nTheta+1)]
fixed = list(map(str, fixed))
print("# group constraint[xyz]", " ".join(fixed))
print("# group tracking", index(0, nL // 2), index(nTheta, nL // 2))
