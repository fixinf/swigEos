import tqdm

shell = get_ipython().__class__.__name__
if shell == 'ZMQInteractiveShell':
    tqdm = tqdm.tqdm_notebook
else:
    tqdm = tqdm.tqdm