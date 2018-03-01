import matplotlib as mpl
import numpy as np
import Label

def plotMRConstraints(obj, alpha=1):
    ax = Label.get_axes(obj)
    # klahnMR1 = np.loadtxt('GrabbedFigures/KlahnMR/KlahnMR1.csv')
    lat_R = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/Lattimer2014/Lattimer2014Right.csv',
        skiprows=1,
        delimiter=','
    )

    lat_L = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/Lattimer2014/Lattimer2014Left.csv',
        skiprows=1,
        delimiter=','
    )

    lat_C = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/Lattimer2014/Lattimer2014Contour.csv',
        skiprows=1,
        delimiter=','
    )

    ax.fill(*lat_C.transpose(), color='#FF9339', fill=0, 
            hatch='\\\\', linewidth=2, zorder=-5, alpha=alpha)
    # ax.plot(*lat_L.tr, anspose())

    ax.fill_between(np.linspace(*ax.get_xlim()), 1.97, 2.05,
             hatch=r'//', zorder=-4, facecolor='None', edgecolor='black',
             alpha=alpha)

    Ozel_LL = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/Ozel16/Left_light.csv',
        skiprows=1,
        delimiter=',')

    ax.fill(*Ozel_LL.transpose(), color='#C9F3F3', fill=1, 
            zorder=-6, alpha=alpha)

    
    Sul68 = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/Suleimanov2016/68p.dat',
        skiprows=0,
        delimiter=' ')

    ax.fill(*Sul68.transpose(), color='#EA83BD', fill=1, 
        zorder=-6, alpha=alpha)

    Bog2s = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/BogdanovMR/Bogdanov_MR2_2s.csv',
        skiprows=1,
        delimiter=' ')

    ax.fill(*Bog2s.transpose(), color='#FFFF70', fill=1, 
        zorder=-17, alpha=alpha)

    KlCaus = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/KlahnMR/KlahnMRCaus.csv',
        skiprows=1,
        delimiter=',')
    
    ax.fill_between(KlCaus[:, 0], KlCaus[:, 1], 2.49,
            zorder=-1, facecolor='white', edgecolor='black',
            alpha=1., lw=0)
    ax.plot(KlCaus[:, 0], KlCaus[:, 1], c='black', ls='-', lw=1)


    Kl1 = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/KlahnMR/KlahnMR1.csv',
        skiprows=1,
        delimiter=',')

    ax.fill_between(Kl1[:, 0], Kl1[:, 1], 1.65077,
            zorder=-12, facecolor='#CAFF70', edgecolor='#099509',
            alpha=alpha, lw=0)#, hatch='XX')

    # ax.text(5.5, 0.9, '4U 0614+09', color='#099509', fontsize=19, backgroundcolor='white')

    Kl2 = np.loadtxt(
        '/home/const/Dropbox/GrabbedFigures/KlahnMR/KlahnMR2.csv',
        skiprows=1,
        delimiter=',')

    ax.fill_between(Kl2[:, 0], Kl2[:, 1], 2.49,
            zorder=-150, facecolor='#C9C9C9', edgecolor='#D3C6C6',
            alpha=alpha, lw=0)#, hatch='\\\\')

    ax.text(7.5, 1.3, '(a)', color='black',
     fontsize=17, backgroundcolor='None')

    ax.text(9.5, 1.8, '(b)', color='black',
    fontsize=17, backgroundcolor='None')

    ax.text(12.5, 1.7, '(c)', color='black',
     fontsize=17, backgroundcolor='None')

    ax.text(13.5, 1.15, '(d)', color='black',
     fontsize=17, backgroundcolor='None')

    ax.text(16, 0.8, '(e)', color='black',
     fontsize=17, backgroundcolor='None')

    ax.text(12, 0.7, '(f)', color='black',
     fontsize=17, backgroundcolor='None')
    # ax.plot(KlCaus[:, 0], KlCaus[:, 1], c='black', ls='-', lw=2)
    # ax.xaxis.set_zorder(6)
    # ax.yaxis.set_zorder(6)
    

def plotMNConstraints(ax):
    pass

