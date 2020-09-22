"""
plotter.py
"""

import numpy as np

def plotter(self, fin, max=1, sym=1):
    """
    Plot function on psgrid
    """

    #Axial connection
    fout = np.squeeze(self.square(fin))
    full = np.hstack( (sym * np.flip(fout, axis=1), fout))
    full = np.vstack( (np.flip(full[0,:]), full, np.flip(full[-1,:])))

    Z = np.squeeze(self.square( self.a * np.cosh(self.Xr) * np.cos(self.Xa)))
    Zf = np.hstack((np.flip(Z, axis=1), Z))
    Zf = np.vstack( (np.flip(Zf[0,:]), Zf, np.flip(Zf[-1,:]) )  )

    X = np.squeeze(self.square( self.a * np.sinh(self.Xr) * np.sin(self.Xa) ))
    Xf = np.hstack( (-np.flip(X, axis=1),X) )
    Xf = np.vstack( (np.flip(Xf[0,:]), Xf, np.flip(Xf[-1,:])) )


    # fig = go.Figure(data=[go.Surface(x=x,y=z,z=full)])

    # fig.update_traces(contours_z=dict(show=True, usecolormap=True,
    #                               highlightcolor="limegreen", project_z=True))

    # fig.update_layout({"template" : "plotly_white"})
    # fig.show()

    return full, Zf, Xf
