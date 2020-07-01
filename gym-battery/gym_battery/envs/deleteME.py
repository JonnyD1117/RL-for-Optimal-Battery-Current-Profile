

from oct2py import octave
from matplotlib import pyplot as plt
import os
import sys

octave.addpath(r"/home/indy/Desktop/PhD/Lin Research/RL Battery Profiler/Github_Repo/RL-for-Optimal-Battery-Current-Profile/gym-battery/gym_battery/envs/")
out1, out2 = octave.model(nout=2)

out1 = out1[0]
out2 = out2[0]


plt.figure(1)
plt.plot(out1, out2)
plt.show()

