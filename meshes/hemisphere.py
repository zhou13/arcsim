import math

nPhi = 240
nTheta = 160

R = 10.

dPhi = 2 * math.pi / nPhi
dTheta = math.pi / 2 / nPhi

for i in range(nPhi):
    for j in range(nTheta):
        phi = dPhi * i
        theta = dTheta * j
        print("v",
              R * math.cos(phi) * math.cos(theta),
              R * math.cos(phi) * math.cos(theta),
              R * math.sin(theta))


def index(i, j):
    return i * (nL + 1) + j


for i in range(nTheta):
    for j in range(nL):
        print("f", index(i, j) + 1, index(i, j+1) + 1, index(i+1, j+1) + 1, index(i+1, j) + 1)

fixed = [index(i, 0) for i in range(nTheta+1)] + [index(i, nL) for i in range(nTheta+1)]
fixed = list(map(str, fixed))
print("# fixed", ", ".join(fixed))
print("# measure", index(0, nL // 2), index(nTheta, nL // 2))
print("# boundary", index(nTheta // 2, 0))
