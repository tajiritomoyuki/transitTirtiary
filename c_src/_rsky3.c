//python関連のモジュールをインポート
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"

//計算関連のモジュールをインポート？
#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

//並列処理用？？
#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

//パイを定義
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

//min関数とmax関数を予め定義？
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

//呼び出せる関数を予め定義？？
static PyObject *_rsky(PyObject *self, PyObject *args);
//static PyObject *_getf(PyObject *self, PyObject *args);

//calculates the eccentric anomaly (see Seager Exoplanets book:  Murray & Correia eqn. 5 -- see section 3)
inline double getE(double M, double e)
{
	double E = M, eps = 1.0e-7;
	double fe, fs;

	// modification from LK 05/07/2017:
	// add fmod to ensure convergence for diabolical inputs (following Eastman et al. 2013; Section 3.1)
	while(fmod(fabs(E - e*sin(E) - M), 2.*M_PI) > eps)
	{
		fe = fmod(E - e*sin(E) - M, 2.*M_PI);
		fs = fmod(1 - e*cos(E), 2.*M_PI);
		E = E - fe/fs;
	}
	return E;
}

static PyObject *_rsky_or_f(PyObject *self, PyObject *args, int f_only)
{
  //変数を定義
  double ecc1, ecc2, inc1, inc2, a1, a2, omega1, omega2, per1, per2, tc1, tc2, BIGD = 100.;
  int transittype, nthreads;
  //これが不明
  npy_intp dims[3];
  //配列を定義？？
  PyArrayObject *ts, *ds12, *ds13, *ds23;
  //pythonからの引数をこちらの変数に代入
  if(!PyArg_ParseTuple(args,"Oddddddii", &ts, &tc1, &per1, &a1, &inc1, &ecc1, &omega1, &tc2, &per2, &a2, &inc2, &ecc2, &omega2, &transittype, &nthreads)) return NULL;
  //dimsにts(時間)を代入
  dims[0] = PyArray_DIMS(ts)[0];
  //1次元dimsと同じ長さ、tsと同じタイプの配列を作成
	ds12 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));
  ds13 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));
  ds23 = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ts));
  //時間の配列を作成
  double *t_array = PyArray_DATA(ts);
  //出力用の配列を作成
	double *array12 = PyArray_DATA(ds12);
  double *array13 = PyArray_DATA(ds13);
  double *array23 = PyArray_DATA(ds23);
  //変数を定義
  //mean motion
  const double n1 = 2.*M_PI/per1;
  const double n2 = 2.*M_PI/per2;
  //小さい値
  const double eps = 1.0e-7;
  //並列計算用。。
	#if defined (_OPENMP) && !defined(_OPENACC)
	omp_set_num_threads(nthreads);	//specifies number of threads (if OpenMP is supported)
	#endif
  //一旦削除
  //#if defined (_OPENACC)
	//#pragma acc parallel loop copyin(t_array[:dims[0]]) copyout(output_array[:dims[0]])
	#if defined (_OPENMP)
	#pragma omp parallel for
	#endif

  //実際に各tごとにrskyを計算する
  for(int i = 0; i < dims[0]; i++)
  {
    //時間
    double t = t_array[i];
    //primaryのtransitが起こる時のtpを求める
    //innerの軌道
		double f1 = M_PI/2. - omega1;
		double E1 = 2.*atan(sqrt((1. - ecc1)/(1. + ecc1))*tan(f1/2.));
		double M1 = E1 - ecc1*sin(E1);
		double tp1 = tc1 - per1*M1/2./M_PI;
    //outerの軌道
		double f2 = M_PI/2. - omega2;
		double E2 = 2.*atan(sqrt((1. - ecc2)/(1. + ecc2))*tan(f2/2.));
		double M2 = E2 - ecc2*sin(E2);
		double tp2 = tc2 - per2*M2/2./M_PI;
    //時刻tのtrue anomalyを計算
    //innerのtrue anomaly
    if(ecc1 < 1.0e-5)
		{
			f1 = ((t - tp1)/per1 - (int)((t - tp1)/per1))*2.*M_PI;
		}
    else
		{
			M1 = n1*(t - tp1);
			E1 = getE(M1, ecc1);
			f1 = 2.*atan(sqrt((1.+ecc1)/(1.-ecc1))*tan(E1/2.));
		}
    //outerのture anomaly
    if(ecc2 < 1.0e-5)
		{
			f2 = ((t - tp2)/per2 - (int)((t - tp2)/per2))*2.*M_PI;
		}
    else
		{
			M2 = n2*(t - tp2);
			E2 = getE(M2, ecc2);
			f2 = 2.*atan(sqrt((1.+ecc2)/(1.-ecc2))*tan(E2/2.));
		}
    //変数を追加
    double d12, d13, d23;
    double x2, x3, y2, y3;
    //場合分けせずにいったんすべて計算する
    //1体目と2体目のrskyを計算
    d12 = a1 * (1.0 - ecc1*cos(E1)) * sqrt(1.0-sin(omega1+f1)*sin(omega1+f1)*sin(inc1)*sin(inc1));
    //1体目と3体目のrskyを計算
    d13 = a1 * (1.0 - ecc1*cos(E2)) * sqrt(1.0-sin(omega2+f2)*sin(omega2+f2)*sin(inc2)*sin(inc2));
    //2体目と3体目のxypositionを計算
    x2 = (cos(omega1) * cos(omega1 + f1) + sin(omega1) * sin(omega1 + f1) * cos(inc1)) * a1 * (1.0 - ecc1*cos(E1));
    y2 = (sin(omega1) * cos(omega1 + f1) + cos(omega1) * sin(omega1 + f1) * cos(inc1)) * a1 * (1.0 - ecc1*cos(E1));
    x3 = (cos(omega2) * cos(omega2 + f2) + sin(omega2) * sin(omega2 + f2) * cos(inc2)) * a2 * (1.0 - ecc2*cos(E2));
    y3 = (sin(omega2) * cos(omega2 + f2) + cos(omega2) * sin(omega2 + f2) * cos(inc2)) * a2 * (1.0 - ecc2*cos(E2));
    //2体目と3体目のrskyを計算
    d23 = sqrt((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3));
    array12[i] = d12;
    array13[i] = d13;
    array23[i] = d23;
    return Py_BuildValue("d", ds12, ds13, ds23);
  }
}

static PyObject *_rsky3(PyObject *self, PyObject *args)
{
	return _rsky_or_f(self, args, 0);
}

static char _rsky3_doc[] = "test";

static PyMethodDef _rsky3_methods[] = {{"_rsky3", _rsky3,METH_VARARGS,_rsky3_doc}};


#if PY_MAJOR_VERSION >= 3
	static struct PyModuleDef _rsky3_module = {
		PyModuleDef_HEAD_INIT,
		"_rsky3",
		_rsky3_doc,
		-1,
		_rsky3_methods
	};

	PyMODINIT_FUNC
	PyInit__rsky(void)
	{
		PyObject* module = PyModule_Create(&_rsky3_module);
		if(!module)
		{
			return NULL;
		}
		import_array();
		return module;
	}
#else
	void init_rsky3(void)
	{
	  Py_InitModule("_rsky3", _rsky3_methods);
	  import_array();
	}
#endif
