"""
The tsnet.utils.calc_parabola_vertex contains function to
calculate the parameters of a parabola based on the
coordinated of three points on the curve.

"""
def calc_parabola_vertex(points):
    """Adapted and modifed to get the unknowns for defining a parabola

    Parameters
    ----------
    points : list
        Three points on the pump characteristics curve.
    """
    [(x1,y1),(x2,y2),(x3,y3)] = points
    denom = (x1-x2) * (x1-x3) * (x2-x3)
    A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom
    B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom
    C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom
    return A,B,C