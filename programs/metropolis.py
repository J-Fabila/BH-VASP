import glob
import os
import math
import random
import re
import sys
import datetime
from subprocess import*

E0 = float(sys.argv[1])
En = float(sys.argv[2])
temperature_K = float(sys.argv[3])

kBT = float(0.00008617 * temperature_K)


if math.exp((E0 - En)/kBT) > random.uniform(0,1):
                print("true")
else:
                print("false")
#print E0 En kBT


