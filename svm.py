from sklearn import svm
import pandas as pd
import numpy as np
import os, threading

def read_mfcc(filename):
  frame = pd.read_csv(filename, sep=' ')
  flat = frame[0:40].as_matrix().flatten()
  return np.pad(flat, (0,14*40 - flat.shape[0]), 'constant', constant_values=(0))

def read_mfccs(dirname):
  res = []
  for filename in os.listdir(dirname):
    if filename.endswith('.mfcc'):
      res.append(read_mfcc(dirname + '/' + filename))
  return res

# hello = read_mfccs('hello')
# other = read_mfccs('other')
# 
# X = np.array(hello + other)
# y = np.array([ 1 for h in hello ] + [ 0 for o in other ])
# 
# clf = svm.LinearSVC()
# 
# clf.fit(X, y)

import subprocess
import plot

def realtime_mfcc():
  proc = subprocess.Popen(['./realtime_mfcc.sh'],stdout=subprocess.PIPE)
  while True:
    line = [ float(s) for s in proc.stdout.readline().split(b' ') ]
    if line != '':
      plot.push(line)
    else:
      break

thd = threading.Thread(target=realtime_mfcc, daemon=True)
thd.start()
plot.render()
