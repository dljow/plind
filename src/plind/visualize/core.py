import numpy as np
import matplotlib.pyplot as plt
# from ..plmodel import PLModel

def visualize_descent(*args, **kwargs):
    pass

# def animate_1d(i):
#     if i == 0:
#         plot_1d(self, step = i, with_background = True, with_contour = True, with_thresh = True)
#     else:
#         plot_1d(self, step = i, with_background = False, with_contour = False, with_thresh = False)


def plot_1d(self, step = -1, with_background = True, with_contour = True, with_thresh = True):

    contour = self.trajectory[step]
    x = contour.points.real.squeeze()
    y = contour.points.imag.squeeze()

    rmax = np.max(np.sqrt(x**2+y**2))
    max_res = np.min((x.size,100))
    axis = np.linspace(-rmax,rmax,max_res)
    XX, YY = np.meshgrid(axis, axis)

    simp = contour.simplices
    simp_lines_x = np.insert(x[simp], 2, np.nan, axis=1)
    simp_lines_y = np.insert(y[simp], 2, np.nan, axis=1)

    # plt.subplot(131)
    if with_background == True:
        vmed = np.median(self.expfun(XX+1j*YY,*self.expargs).real)
        vmax = vmed+np.median(np.abs(self.expfun(XX+1j*YY,*self.expargs).real-vmed))*4
        vmin = vmed-np.median(np.abs(self.expfun(XX+1j*YY,*self.expargs).real-vmed))*4
        plt.imshow(np.flipud(self.expfun(XX+1j*YY,*self.expargs).real),vmin=vmin,vmax=vmax,aspect='auto',interpolation='bicubic',extent=(XX.min(),XX.max(),YY.min(),YY.max()))
    if with_contour == True:
        imed = np.median(self.expfun(XX+1j*YY,*self.expargs).imag)
        imax = imed+np.median(np.abs(self.expfun(XX+1j*YY,*self.expargs).imag-imed))*4
        imin = imed-np.median(np.abs(self.expfun(XX+1j*YY,*self.expargs).imag-imed))*4
        levels = np.linspace(imin,imax,10)
        print(levels)
        plt.contour(XX,YY,self.expfun(XX+1j*YY,*self.expargs).imag,levels=levels,cmap=plt.cm.magma,linewidths=0.7)
    if with_thresh == True:
        plt.contour(XX,YY,self.expfun(XX+1j*YY,*self.expargs).real, levels=[self.thresh], colors='black')
        plt.contourf(XX,YY,self.expfun(XX+1j*YY,*self.expargs).real, colors='black',levels=[-np.inf,self.thresh],alpha=0.0,rasterized=True,hatches=['//',None])

    simp_lines = plt.plot(simp_lines_x.ravel(), simp_lines_y.ravel(),color='tab:red')
    plt.plot(x,y,'.',color='r',markersize=1.5)
    plt.xlabel('Re(x)')
    plt.ylabel('Im(x)')

    plt.tight_layout()
    return
