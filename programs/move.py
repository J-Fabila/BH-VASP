import random
import sys
step_width = float(sys.argv[1])

move = 2 * step_width * random.uniform(-1, 1)

print(move)
