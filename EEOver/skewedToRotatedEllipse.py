from matplotlib import pyplot as plt
from numpy import arctan as atan
from numpy import tan as tan, cos
from numpy import sqrt
from matplotlib.patches import Ellipse

def skewedToRotatedEllipse(width, height, skew):
    c = skew;
    hs = height;
    ws = width;
    if c == 0:
        phi = 0.0
        w = width/2
        h = height/2
    else:
        phi = atan((c**2*hs**2 - hs**2 + ws**2 + sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))/(2*c*hs**2));
        w = sqrt(hs**2*ws**2*(4*c**2*hs**4 + (c**2*hs**2 - hs**2 + ws**2 + sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2)/(4*c**2*hs**6 + 2*c**2*hs**4*(4*c**2*hs**4 + (c**2*hs**2 - hs**2 + ws**2 + sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2)*(c**2*hs**2 - hs**2 + ws**2 + sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4))/(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + c**2*hs**2*sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4) + hs**4 - 2*hs**2*ws**2 - hs**2*sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4) + ws**4 + ws**2*sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4)) + c**2*hs**2*(c**2*hs**2 - hs**2 + ws**2 + sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2 + ws**2*(c**2*hs**2 - hs**2 + ws**2 + sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2))/2;
        h =  sqrt(hs**2*ws**2*(4*c**2*hs**4 + (c**2*hs**2 - hs**2 + ws**2 - sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2)/(4*c**2*hs**6 + 2*c**2*hs**4*(4*c**2*hs**4 + (c**2*hs**2 - hs**2 + ws**2 - sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2)*(c**2*hs**2 - hs**2 + ws**2 - sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4))/(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 - c**2*hs**2*sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4) + hs**4 - 2*hs**2*ws**2 + hs**2*sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4) + ws**4 - ws**2*sqrt(c**4*hs**4 + 2*c**2*hs**4 + 2*c**2*hs**2*ws**2 + hs**4 - 2*hs**2*ws**2 + ws**4)) + c**2*hs**2*(c**2*hs**2 - hs**2 + ws**2 - sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2 + ws**2*(c**2*hs**2 - hs**2 + ws**2 - sqrt((c**2*hs**2 + hs**2 - 2*hs*ws + ws**2)*(c**2*hs**2 + hs**2 + 2*hs*ws + ws**2)))**2))/2;
    return phi, w, h


def rotatedToSkewedEllipse(ws, hs, phi):
    if phi == 0:
        skew = 0.0
    else:
        skew = (hs**2*tan(phi)**2 - hs**2 - sqrt(hs**4/cos(phi)**4 - 4*hs**2*ws**2*tan(phi)**2))/(2*hs**2*tan(phi))
    return skew

if __name__ == "__main__":
    import numpy as np
    X1 = 0.0
    Y1 = 0.0

    fig, axes = plt.subplots(1,5, sharex=True,sharey=True, figsize=(20,5))
    for height in np.linspace(0.1,2.0,5):
        width = height;
        for i, skew in enumerate(np.linspace(-1.0,1.0,5)):
            (phi,w,h) = skewedToRotatedEllipse(width, height, skew)
            ellipse = Ellipse(xy=(X1,Y1), width=2.*w, height=2.*h, angle=phi*180./np.pi, edgecolor='r', fc='None')
            axes[i].add_patch(ellipse)
            axes[i].set_title('skew %.3f' % skew)
            axes[i].set_xlim([-5.0,5.0])
            axes[i].set_ylim([-5.0,5.0])
            axes[i].plot(-(height/2)*skew,(height/2),'x')
            axes[i].axis('equal')

    plt.show()

