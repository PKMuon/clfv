import multiprocessing
import os

hostname = os.getenv('HOSTNAME')
nproc = multiprocessing.cpu_count()
if 'phy.pku.edu.cn' in hostname: nproc //= 2
pool = multiprocessing.Pool(nproc)
