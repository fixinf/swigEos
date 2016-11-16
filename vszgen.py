import veusz.embed as vsz
from os.path import join, isfile
import numpy as np
import uuid

VSZ_TEMPLATE_DIR = '/home/const/Dropbox/vsz_templates/'
_templates = {'mr' : 'mr.vsz',
              'mn' : 'mn.vsz',
              'pn' : 'pn.vsz',
              'pods' : 'pods.vsz'}

VSZ_EXPORT_DIR = '/home/const/Dropbox/vsz_export/'

templates = {t : join(VSZ_TEMPLATE_DIR, _templates[t]) for t in _templates}

class vszSimplePlot(vsz.Embedded):
    def __init__(self, hidden=1, template=None, filename=None):
        super().__init__(hidden=hidden)
        self.lines = []
        if template:
            if template in templates:
                self.Load(templates[template])
            else:
                print('No template "' + template + '" exists; aborting')
                raise
        else:
            self.page = self.Root.Add('page')
            self.graph = self.page.Add('graph')

        self.colors = ['red', 'blue', 'green', 'magenta', 'cyan']
        self.styles = ['solid', 'dashed', 'dot1', 'dash-dot', 'dash-dot-dot']
        self.drawn_labels = []
        self.filename = filename

    def ls_Main(self, l, color=None, style=None):
        l.MarkerLine.hide.val = True
        l.MarkerFill.hide.val = True
        l.PlotLine.width.val = '3pt'
        if color:
            l.PlotLine.color.val = color
        if style:
            l.PlotLine.style.val = style

    def ls_Thin(self, l, color=None, style=None):
        self.ls_Main(l, color=color, style=style)
        l.PlotLine.width.val = '1.5pt'

    def drawLine(self, x, y, label, num_lines=5, colors=None, cycle=1, name=0):
        if label in self.drawn_labels:
            raise ValueError('Line label must be unique!')

        if name:
            xname = label + '_x'
            yname = label + '_y'
            l = self.graph.Add('xy', name=label, index=0)
        else:
            xname = str(uuid.uuid4())
            yname = str(uuid.uuid4())
            l = self.graph.Add('xy', index=0)

        self.SetData(xname, x)
        self.SetData(yname, y)

        l.xData.val = xname
        l.yData.val = yname
        l.key.val = label
        if cycle:
            self.lines.append(l)
            self.drawn_labels.append(label)

        self.ls_Main(l, color=self.colors[:num_lines][::-1][len(self.lines) - 1],
                        style=self.styles[:num_lines][::-1][len(self.lines) - 1])

        return l

    def SaveCheck(self, filename, force=0):
        path = join(VSZ_EXPORT_DIR, filename)
        if isfile(path) and not force:
            print('File ' + path + ' already exists! Not written')
            return
        else:
            # print('File ' + path + ' doesn\'t exist; writing')
            self.Save(path, mode='vsz')


class MassPlot(vszSimplePlot):
    def __init__(self, models, labels, template=None, hidden=1, filename=None):
        super(MassPlot, self).__init__(hidden=hidden, template=template,
                                    filename=filename)
        self.graph = self.Root['page1']['graph1']
        self.models = models
        self.labels = labels
        self.x_row = None
        self.y_row = None

    def drawLines(self):
        #Check what to plot
        if self.x_row is None or self.y_row is None:
            raise NotImplementedError('Proper x_row and y_row \
            should be defined in the constructor!')

        #Plot in inverse order for conventional line overlap
        for model, label in list(zip(self.models, self.labels))[::-1]:
            try:
                data = np.loadtxt(join(model.foldername,
                                   model.filenames['mass_crust']+'_linear'),
                                   skiprows=1)
            except FileNotFoundError:
                raise FileNotFoundError('No data found for model ' + label)
                return

            #Index of a maximum mass -- fixed row #2
            imax = np.argmax(data[:, 1])

            self.drawLine(data[:imax, self.x_row], data[:imax, self.y_row],
             label, num_lines=len(self.models))

            self.drawMaxEllipse(data[:, x_row], data[:, y_row],
                self.colors[:len(self.models)][::-1][len(self.lines) - 1])

    def drawMaxEllipse(self, x, y, color):
        imax = np.argmax(y)
        xmax, ymax = x[imax], y[imax]
        e = self.graph.Add('ellipse', index=0)
        e.width.val = 0.03
        e.height.val = 0.03
        e.Border.hide.val = True
        e.Fill.hide.val = False
        e.Fill.color.val = color
        e.positioning.val = 'axes'
        e.xPos.val = xmax
        e.yPos.val = ymax

class MRPlot2(MassPlot):
    def __init__(self, models, labels, hidden=1, filename=None):
        super(MRPlot2, self).__init__(models, labels, hidden=hidden,
                        template='mr', filename=filename)

        self.y_row = 1
        self.x_row = 2
        self.drawLines()
        if self.filename:
            self.SaveCheck(self.filename)

class MRPlot(vszSimplePlot):
    """docstring for MRPlot."""
    def __init__(self, models, labels, hidden=1, filename=None):
        super(MRPlot, self).__init__(hidden=hidden, template='mr',
                                    filename=filename)
        self.graph = self.Root['page1']['graph1']
        for model, label in zip(models[::-1], labels[::-1]):
            try:
                data = np.loadtxt(join(model.foldername,
                                   model.filenames['mass_crust']+'_linear'),
                                   skiprows=1)
            except FileNotFoundError:
                print('no data')
                return

            imax = np.argmax(data[:, 1])

            self.drawLine(data[:imax, 2], data[:imax, 1], label,
                num_lines=len(models))
            self.drawMaxEllipse(data[:, 2], data[:, 1],
                self.colors[:len(models)][::-1][len(self.lines) - 1])

        if self.filename:
            self.SaveCheck(self.filename)

    def drawMaxEllipse(self, x, y, color):
        imax = np.argmax(y)
        xmax, ymax = x[imax], y[imax]
        e = self.graph.Add('ellipse', index=0)
        e.width.val = 0.03
        e.height.val = 0.03
        e.Border.hide.val = True
        e.Fill.hide.val = False
        e.Fill.color.val = color
        e.positioning.val = 'axes'
        e.xPos.val = xmax
        e.yPos.val = ymax

class MNPlot(vszSimplePlot):
    """docstring for MRPlot."""
    def __init__(self, models, labels, hidden=1, filename=None):
        super(MNPlot, self).__init__(hidden=hidden, template='mn',
                                    filename=filename)
        self.graph = self.Root['page1']['graph1']
        for model, label in zip(models[::-1], labels[::-1]):
            try:
                data = np.loadtxt(join(model.foldername,
                                   model.filenames['mass_crust']+'_linear'),
                                   skiprows=1)
            except FileNotFoundError:
                print('no data')
                return

            imax = np.argmax(data[:, 1])
            self.drawLine(data[:imax, 0], data[:imax, 1], label,
                num_lines=len(models))
            self.drawMaxEllipse(data[:, 0], data[:, 1],
                self.colors[:len(models)][::-1][len(self.lines) - 1])

        if self.filename:
            self.Save(join(VSZ_EXPORT_DIR, self.filename), mode='vsz')

    def drawMaxEllipse(self, x, y, color):
        return MRPlot.drawMaxEllipse(self, x, y, color)


class RNPlot(vszSimplePlot):
    def __init__(self, models, labels, hidden=1, filename=None):
        super(RNPlot, self).__init__(hidden=hidden, template='rn',
                                    filename=filename)
        self.graph = self.Root['page1']['graph1']
        for model, label in zip(models[::-1], labels[::-1]):
            try:
                data = np.loadtxt(join(model.foldername,
                                   model.filenames['mass_crust']+'_linear'),
                                   skiprows=1)
            except FileNotFoundError:
                print('no data')
                return

            imax = np.argmax(data[:, 1])
            self.drawLine(data[:imax, 0], data[:imax, 1], label,
                num_lines=len(models))
            self.drawMaxEllipse(data[:, 0], data[:, 1],
                self.colors[:len(models)][::-1][len(self.lines) - 1])

        if self.filename:
            self.Save(join(VSZ_EXPORT_DIR, self.filename), mode='vsz')


class PNPlot(vszSimplePlot):
    def __init__(self, models, invs, labels, hidden=1, filename=None,
        xlim=None, ylim=None):
        super(PNPlot, self).__init__(hidden=hidden, template='pn',
                                    filename=filename)
        self.graph = self.Root['page1']['graph1']
        for model, label, inv in list(zip(models, labels, invs))[::-1]:
            try:
                if inv:
                    model.loadEosInv()
                    self.drawLine(model.nrange_inv/model.n0,
                                model._P_inv, label,
                                num_lines=len(models), cycle=1)
                    #
                    if model.needsMaxwInv():
                        model.processMaxwInv(shift=200)
                        lmaxw = self.drawLine(model.nrange_maxw/model.n0,
                        model._P_maxw, '', cycle=0)
                        self.ls_Thin(lmaxw,
                        color=self.lines[-1].PlotLine.color.val,
                        style=self.lines[-1].PlotLine.style.val)



                else:
                    model.loadEos()
                    self.drawLine(model.nrange/model.n0,
                                model._P, label,
                        num_lines=len(models))
            except FileNotFoundError:
                print('No data for model ' + label)
                # return
