__author__ = 'const'
import multiprocessing as mp

def f(x):
    print("Hello, world!", mp.current_process())
    return 1

pool = mp.Pool(4)


pool.map(f, data)