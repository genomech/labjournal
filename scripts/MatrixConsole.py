import os

shape = os.popen('stty size', 'r').read().split()
shape = [int(x) for x in shape]
screen = [[0] * shape[1]] * shape[0]

print(screen)
