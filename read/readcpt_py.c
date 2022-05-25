#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "readcpt.c"

static int setpixeldict(PyObject *pixeldict, struct cpt_pixel *ppixel);

/*  Main fn src  */
static PyObject *cpt_readall_py(PyObject *self, PyObject *args)
{
	char *fname = NULL;
	uint8_t  nparam, iparam, ipoint, ivicinity;
	uint32_t nptx, iptx;
	struct cpt_ptx ptx;
	struct cpt_pt *ppt;
	struct cpt_px *ppx;
	struct cpt_point *ppoint;
	struct cpt_pixel *ppixel;
	PyObject *retlist,
	         *ptxdict,
	         *ptdict,
	         *pxdict,
	         *pointlist,
	         *pointdict,
	         *datalist,     /*  data of site parameters  */
	         *pixeldict,
	         *vicilist;
	
	/*  Wrap pystring to char*  */
	if(!PyArg_ParseTuple(args, "s", &fname))
		return NULL;
	
	/*  Original C result  */
	if ((cpt_readall(fname, &ptx, &nptx, &nparam)))
		return NULL;
	
	/*  Wrap C result to python  */
	retlist = PyList_New(nptx);
	for (iptx = 0; iptx < nptx; ++iptx) {
		ptxdict = PyDict_New();
		
		/*  Pt  */
		ppt = ptx.pt+iptx;
		ptdict = PyDict_New();
		pointlist = PyList_New(ppt->nt);
		for (ipoint = 0; ipoint < ppt->nt; ++ipoint) {
			ppoint = ppt->points+ipoint;
			
			pointdict = PyDict_New();
			datalist  = PyList_New(nparam);
			for (iparam = 0; iparam < nparam; ++iparam) {
				PyList_SetItem(datalist, iparam,
				               PyFloat_FromDouble(ppoint->params[iparam]));
			}
			PyDict_SetItemString(pointdict, "data", datalist);
			PyDict_SetItemString(pointdict, "time", PyLong_FromLongLong(ppoint->seconds));
			
			PyList_SetItem(pointlist, ipoint, pointdict);
		}
		PyDict_SetItemString(ptdict, "point", pointlist);
		PyDict_SetItemString(ptdict, "lon", PyFloat_FromDouble(ppt->lon));
		PyDict_SetItemString(ptdict, "lat", PyFloat_FromDouble(ppt->lat));
		PyDict_SetItemString(ptdict, "alt", PyLong_FromLong(ppt->alt));
		
		/*  Px  */
		ppx = ptx.px+iptx;
		pxdict = PyDict_New();
		PyDict_SetItemString(pxdict, "time", PyLong_FromLongLong(ppx->seconds));
		
		ppixel = ppx->centrepixel;
		pixeldict = PyDict_New();
		setpixeldict(pixeldict, ppixel);
		PyDict_SetItemString(pxdict, "center", pixeldict);
		
		vicilist = PyList_New(ppx->nvicinity);
		for (ivicinity = 0; ivicinity < ppx->nvicinity; ++ivicinity) {
			ppixel = ppx->vicinity+ivicinity;
			pixeldict = PyDict_New();
			setpixeldict(pixeldict, ppixel);
			
			PyList_SetItem(vicilist, ivicinity, pixeldict);
		}
		PyDict_SetItemString(pxdict, "vicinity", vicilist);
		
		PyDict_SetItemString(ptxdict, "pt", ptdict);
		PyDict_SetItemString(ptxdict, "px", pxdict);
		PyList_SetItem(retlist, iptx, ptxdict);
	}
	
	/*  Free original C result  */
	cpt_freeptall(&ptx.pt, nptx);
	cpt_freepxall(&ptx.px, nptx);
	
	/*  Return to python  */
	return retlist;
}

static int setpixeldict(PyObject *pixeldict, struct cpt_pixel *ppixel)
{
	uint8_t  ilayer;
	PyObject *banddict, *channeldict, *satlist;
	struct cpt_channel *pchannel;
	
	banddict = PyDict_New();
	for (uint8_t ichannel = 0; ichannel < ppixel->nchannel; ++ichannel) {
		pchannel = ppixel->channels+ichannel;
		channeldict = PyDict_New();
		
		/*  I  */
		satlist = PyList_New(ppixel->nlayer);
		for (ilayer = 0; ilayer < ppixel->nlayer; ++ilayer) {
			PyList_SetItem(satlist, ilayer,
			               PyFloat_FromDouble(pchannel->obs[ilayer]));
		}
		PyDict_SetItemString(channeldict, "I", satlist);
		
		/*  Q and U  */
		if (pchannel->centrewv < 0) {
			satlist = PyList_New(ppixel->nlayer);
			for (ilayer = ppixel->nlayer; ilayer < 2*ppixel->nlayer; ++ilayer) {
				PyList_SetItem(satlist, ilayer-ppixel->nlayer,
					       PyFloat_FromDouble(pchannel->obs[ilayer]));
			}
			PyDict_SetItemString(channeldict, "Q", satlist);
			
			satlist = PyList_New(ppixel->nlayer);
			for (ilayer = 2*ppixel->nlayer; ilayer < 3*ppixel->nlayer; ++ilayer) {
				PyList_SetItem(satlist, ilayer-2*ppixel->nlayer,
					       PyFloat_FromDouble(pchannel->obs[ilayer]));
			}
			PyDict_SetItemString(channeldict, "U", satlist);
		}
		
		/*  sz/vz/sa/va  */
		satlist = PyList_New(ppixel->nlayer);
		for (ilayer = 0; ilayer < ppixel->nlayer; ++ilayer) {
			PyList_SetItem(satlist, ilayer,
			               PyFloat_FromDouble(pchannel->ang[ilayer]));
		}
		PyDict_SetItemString(channeldict, "sza", satlist);
		
		satlist = PyList_New(ppixel->nlayer);
		for (ilayer = ppixel->nlayer; ilayer < 2*ppixel->nlayer; ++ilayer) {
			PyList_SetItem(satlist, ilayer-ppixel->nlayer,
			               PyFloat_FromDouble(pchannel->ang[ilayer]));
		}
		PyDict_SetItemString(channeldict, "vza", satlist);
		
		satlist = PyList_New(ppixel->nlayer);
		for (ilayer = 2*ppixel->nlayer; ilayer < 3*ppixel->nlayer; ++ilayer) {
			PyList_SetItem(satlist, ilayer-2*ppixel->nlayer,
			               PyFloat_FromDouble(pchannel->ang[ilayer]));
		}
		PyDict_SetItemString(channeldict, "saa", satlist);
		
		satlist = PyList_New(ppixel->nlayer);
		for (ilayer = 3*ppixel->nlayer; ilayer < 4*ppixel->nlayer; ++ilayer) {
			PyList_SetItem(satlist, ilayer-3*ppixel->nlayer,
			               PyFloat_FromDouble(pchannel->ang[ilayer]));
		}
		PyDict_SetItemString(channeldict, "vaa", satlist);
		
		PyDict_SetItem(banddict,
		               PyLong_FromLong(Py_ABS(pchannel->centrewv)), channeldict);
	}
	PyDict_SetItemString(pixeldict, "band", banddict);
	
	PyDict_SetItemString(pixeldict, "lon", PyFloat_FromDouble(ppixel->lon));
	PyDict_SetItemString(pixeldict, "lat", PyFloat_FromDouble(ppixel->lat));
	PyDict_SetItemString(pixeldict, "alt", PyLong_FromLong(ppixel->alt));
	PyDict_SetItemString(pixeldict, "flag", PyLong_FromLong(ppixel->mask));
	
	return 0;
}


/*  Register fn to python  */
static PyMethodDef cptreadallpymethod[] = {
	{"load", cpt_readall_py, METH_VARARGS, "Load entire cpt"},
	{NULL, NULL, 0, NULL}
};
static struct PyModuleDef cptreadallpymod = {
	PyModuleDef_HEAD_INIT,
	"pycpt",
	"Python interface for reading cpt file format file",
	-1,
	cptreadallpymethod
};
PyMODINIT_FUNC PyInit_pycpt(void)
{
	return PyModule_Create(&cptreadallpymod);
}

