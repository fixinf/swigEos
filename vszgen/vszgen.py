import veusz

def createMR(filename, models, labels=[], show=False):
    g = veusz.Embedded(hidden=not show)
    page = g.Root.Add(page)
    graph = page.Add('graph')
    for m in Models:
        xy = graph.Add()
