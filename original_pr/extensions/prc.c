#include <math.h>
#include "Python.h"
#include "numpy/arrayobject.h"

/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/

PyObject *prc(PyObject *self, PyObject *args){
	PyArrayObject *array = NULL ;
    PyObject *pylist, *item ;

	//PyArrayObject *dist ;
	double x1,y1,z1,x2,y2,z2,sdist;
	int i,j,k,natoms,npairs,check;
	
	if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array))
		return NULL;

	natoms = array->dimensions[0];
    npairs = (natoms * (natoms - 1))/2 ;
    double dist[npairs] ;

	check=0 ;
    k=0;
	for (i=0; i< natoms-1 ; i++){
		x1=*(double *)(array->data + i*array->strides[0]+0*array->strides[1]) ;
		y1=*(double *)(array->data + i*array->strides[0]+(0+1)*array->strides[1]) ;
		z1=*(double *)(array->data + i*array->strides[0]+(0+2)*array->strides[1]) ;
		for (j=i+1; j<natoms ; j++){
			x2=*(double *)(array->data + j*array->strides[0]+0*array->strides[1]) ;
			y2=*(double *)(array->data + j*array->strides[0]+(0+1)*array->strides[1]) ;
			z2=*(double *)(array->data + j*array->strides[0]+(0+2)*array->strides[1]) ;
			sdist=((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)) ;
			dist[k]=sqrt(sdist) ;
            k++ ;
			}
		}

    pylist = PyList_New(npairs) ;
    if (pylist != NULL){
        for (i=0 ; i<npairs ; i++) {
            item = PyFloat_FromDouble(dist[i]);
            PyList_SET_ITEM(pylist, i, item);
        }
    }

    return pylist ;

}

static PyMethodDef methods[] = {
	{ "prc", prc, METH_VARARGS },
	{ NULL, NULL }
} ;

void initprc(){
	PyObject *m ;
	m = Py_InitModule("prc", methods);
	import_array();
}

