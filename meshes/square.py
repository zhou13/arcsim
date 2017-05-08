import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--nxy', type=str, help="two numbers seperated by a comma")
parser.add_argument('--xy', type=str, help="two numbers seperated by a comma")
parser.add_argument('--fixed', type=str, help="numbers seperated by a comma")
args = parser.parse_args()

nX, nY = map(int, args.nxy.split(','))
X, Y = map(float, args.xy.split(','))
nfix = list(map(int, args.fixed.split(',')))

print(f"# nX = {nX}\n# nY = {nY}\n# X = {X}\n# Y = {Y}\n")

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
        print("f", index(i, j) + 1, index(i+1, j) + 1, index(i+1, j+1) + 1)
        print("f", index(i, j) + 1, index(i+1, j+1) + 1, index(i, j+1) + 1)

fixed = [index(i, j) for i in range(nX+1) for j in nfix]
fixed = list(map(str, fixed))
print("# group constraint_xyz", ", ".join(fixed))

measure = [index(i, nY) for i in range(nX+1)]
measure = list(map(str, measure))
print("# group tracking", ", ".join(measure))
print("# group force", ", ".join(measure))
# print("# group force ", index(nX-1, nY))
