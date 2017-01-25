import math

nX = 9
nY = 6000

X = 2
Y = 450.

dX = X / nX
dY = Y / nY
for i in range(-nX // 2, nX // 2 + 1):
    for j in range(-nY // 2, nY // 2 + 1):
        x = dX * i
        y = dY * j
        print("v", x, y, 0)


def index(i, j):
    return i * (nY + 1) + j


for i in range(nX):
    for j in range(nY):
        print("f", index(i, j) + 1, index(i+1, j) + 1, index(i+1, j+1) + 1, index(i, j+1) + 1)

fixed = [index(i, j) for i in range(nX+1) for j in range(6)]
fixed = list(map(str, fixed))
print("# fixed", ", ".join(fixed))

measure = [index(i, nY) for i in range(nX+1)]
measure = list(map(str, measure))
print("# measure", ", ".join(measure))
