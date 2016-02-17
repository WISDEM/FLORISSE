#include "solvers.h"
#include "program_constants.h"
#include <Python.h>


static PyObject* ellipse_ellipse_overlap(PyObject* self, PyObject* args)
{
    int nroots_gsl = 0;
    int rtn_gsl = 0;

    int choice = 1;
    double x_gsl[4], y_gsl[4];

    double A1;
    double B1;
    double H1;
    double K1;
    double PHI_1;
    double A2;
    double B2;
    double H2;
    double K2;
    double PHI_2;

    double area_gsl;

    if (!PyArg_ParseTuple(args, "dddddddddd", &A1, &B1, &H1, &K1, &PHI_1, &A2, &B2, &H2, &K2, &PHI_2))
        return NULL;

    area_gsl = ellipse_ellipse_overlap_gsl(PHI_1, A1, B1, H1, K1, PHI_2, A2, B2, H2, K2, x_gsl, y_gsl, &nroots_gsl, &rtn_gsl, choice);

    return Py_BuildValue("f", area_gsl);
}


static PyMethodDef eeMethods[] = 
{ 
    { "ellipse_ellipse_overlap", ellipse_ellipse_overlap, METH_VARARGS, "ellipse-ellipse overlap"}, 
    { NULL, NULL, 0, NULL }
};

/* module initialization */
PyMODINIT_FUNC 
initellipseEllipseOverlap(void)
{
   (void) Py_InitModule("ellipseEllipseOverlap", eeMethods);
}