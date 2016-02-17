from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
from ellipseEllipseOverlap import ellipse_ellipse_overlap
import numpy as np
import sympy as sp

# plot a skewed ellipse
ws = sp.Symbol("ws")
hs = sp.Symbol("hs")

c = sp.Symbol("c") # skew
r = sp.Symbol("r")
phi = sp.Symbol("phi")

x = r*sp.cos(phi)
y = r*sp.sin(phi)

skewed_ellipse = ((x+c*y)/(ws/2))**2+(y/(hs/2))**2 - 1
r_s = sp.solve(skewed_ellipse,r)[0]
print 'r_s = ', r_s 
diff_r_s = sp.simplify(sp.diff(r_s,phi))
print 'diff_r_s = ', diff_r_s 
min_max_phi = sp.solve(diff_r_s,phi)
PHI_expr = min_max_phi[1]
w_expr = sp.simplify(-1*r_s.subs([(phi,min_max_phi[1])]))
h_expr = sp.simplify(-1*r_s.subs([(phi,min_max_phi[0])]))

PHI = sp.Symbol("PHI")
s_expr = sp.solve(PHI_expr-PHI,c)

print 'PHI = ', PHI_expr
print 's = ', s_expr
print 'w = ', w_expr
print 'h = ', h_expr

#substitute
ws_v = 1.0
hs_v = 2.0
c_v = 0.1

PHI = PHI_expr.subs([(hs,hs_v), (ws,ws_v), (c,c_v)])
w = w_expr.subs([(hs,hs_v), (ws,ws_v), (c,c_v)])
h = -1*r_s.subs([(phi,min_max_phi[0]), (hs,hs_v), (ws,ws_v), (c,c_v)])

print 'PHI', PHI
print 'w', w
print 'h', h

phi_range = np.linspace(0,2*np.pi,30)
r_v = np.array([r_s.subs([(phi,phi_v), (hs,hs_v), (ws,ws_v), (c,c_v)]).evalf() for phi_v in phi_range])
x_v = r_v*np.cos(phi_range)
y_v = r_v*np.sin(phi_range)
plt.plot(x_v,y_v)

rotated_ellipse = (x*sp.cos(PHI)+y*sp.sin(PHI))**2/w**2 + (x*-sp.sin(PHI)+y*sp.cos(PHI))**2/h**2 - 1
r_s = sp.solve(rotated_ellipse,r)[0]
phi_range = np.linspace(0,2*np.pi,30)
r_v = np.array([r_s.subs(phi,phi_v) for phi_v in phi_range])
x_v = r_v*np.cos(phi_range)
y_v = r_v*np.sin(phi_range)
plt.plot(x_v,y_v,'k--')
plt.axis('equal')
plt.show()

#plt.show(block=False)


# t = sp.atan((h/sp.tan(PHI))/w);
# x_y_max = w*sp.cos(t)*sp.cos(PHI) - h*sp.sin(t)*sp.sin(PHI);
# r = sp.solve(x_y_max-c,r)


# # input parameters
# R = sp.Symbol("R")
# phi = sp.Symbol("phi")

# # parameter to solve for
# w = sp.Symbol("w")

# h = (sp.pi*R*R)/(sp.pi*w) # setting height such that area of ellipse is equal to area of circle with radius R

# # find expression for maximum height
# t = sp.atan((h/sp.tan(phi))/w);
# y_max = h*sp.sin(t)*sp.cos(phi) + w*sp.cos(t)*sp.sin(phi);
# x_y_max = w*sp.cos(t)*sp.cos(phi) - h*sp.sin(t)*sp.sin(phi);

# x_y_max

# # find w such that y_max == R
# X1 = 0.
# Y1 = 0.
# PHI_1 = 45.*np.pi/180.
# Rv = 1.0

# y_max_s = y_max.subs(phi, PHI_1)
# y_max_s = y_max_s.subs(R, Rv)

# x_y_max_s = x_y_max.subs(phi, PHI_1)
# x_y_max_s = x_y_max_s.subs(R, Rv)

# h_s = h.subs(R, Rv)

# w_range = np.linspace(0.2,1.0,10)
# y_max_w = [y_max_s.subs(w,w_value) for w_value in w_range]
# x_y_max_w = [x_y_max_s.subs(w,w_value) for w_value in w_range]
# h_w = [h_s.subs(w,w_value) for w_value in w_range]

# fig, axes = plt.subplots(1,10, sharex=True,sharey=True)
# for i, W1 in enumerate(w_range):
#     H1 = h_w[i];
#     ellipse = Ellipse(xy=(X1,Y1), width=2.*W1, height=2.*H1, angle=PHI_1*180./np.pi, edgecolor='r', fc='None')
#     axes[i].add_patch(ellipse)
#     axes[i].plot(x_y_max_w[i],y_max_w[i],'*')
#     axes[i].axis('equal')
#     axes[i].set_ylim([-10,10])
# plt.show()



# W1 = 1.
# H1 = 2.
# X1 = 0.
# Y1 = 0.
# PHI_1 = 45.*np.pi/180.

# W2 = 1.
# H2 = 1.
# X2 = 0.
# Y2 = 0.
# PHI_2 = 45.*np.pi/180.

# overlap_area = ellipse_ellipse_overlap(W1, H1, X1, Y1, PHI_1, W2, H2, X2, Y2, PHI_2)

# plt.figure()
# ax = plt.gca()
# ellipse = Ellipse(xy=(X1,Y1), width=2.*W1, height=2.*H1, angle=PHI_1*180./np.pi, edgecolor='r', fc='None')
# ax.add_patch(ellipse)
# ellipse = Ellipse(xy=(X2,Y2), width=2.*W2, height=2.*H2, angle=PHI_2*180./np.pi, edgecolor='b', fc='None')
# ax.add_patch(ellipse)
# plt.axis('equal')

# # calculate height of ellipse

# t = np.arctan((H1/np.tan(PHI_1))/W1);
# print t
# y_max = H1*np.sin(t)*np.cos(PHI_1) + W1*np.cos(t)*np.sin(PHI_1);
# x_y_max = W1*np.cos(t)*np.cos(PHI_1) - H1*np.sin(t)*np.sin(PHI_1);
# print y_max

# plt.plot(x_y_max,y_max,'*')

# plt.show()