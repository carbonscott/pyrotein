#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# Load upstream data...
u     = np.load("u.seq.npy")
s     = np.load("s.seq.npy")
vh    = np.load("vh.seq.npy")


top = 20 + 1
np.save("u.seq.trunc.npy",  u[:, :top])
np.save("s.seq.trunc.npy",  s[:top])
np.save("vh.seq.trunc.npy", vh[:top])
