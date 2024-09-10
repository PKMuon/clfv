import multiprocessing
import getpass

user = getpass.getuser()
nproc = multiprocessing.cpu_count()
if 'phy.pku.edu.cn' in user: nproc //= 2
pool = multiprocessing.Pool(nproc)
