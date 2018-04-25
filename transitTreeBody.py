#-*-coding:utf-8-*-
import batman
import numpy as np
import multiprocessing
#from . import openmp
from PyAstronomy import pyasl
from . import _rsky

class TransitTreeBody():
    def __init__(self, params, t, nthreads = 1):
        #パラメータの設定
        #角度はすべてdegreeにしよう
        self.t = t
        self.t0_1 = params.t0_1
        self.per_1 = params.per_1
        self.rp_12 = params.rp_12
        self.a_1 = params.a_1
        self.inc_1 = params.inc_1
        self.ecc_1 = params.ecc_1
        self.w_1 = params.w_1
        self.u_1 = params.u_1
        self.limb_dark_1 = params.limb_dark_1
        self.fp_1 = params.fp_1
        self.t_secondary_1 = params.t_secondary_1
        self.t0_2 = params.t0_2
        self.per_2 = params.per_2
        self.rp_13 = params.rp_13
        self.a_2 = params.a_2
        self.inc_2 = params.inc_2
        self.ecc_2 = params.ecc_2
        self.w_2 = params.w_2
        self.u_2 = params.u_2
        self.limb_dark_2 = params.limb_dark_2
        self.fp_2 = params.fp_2
        self.t_secondary_2 = params.t_secondary_2
        self.rp_23 = self.rp_13 / self.rp_23

    #スレッド
    if nthreads==None or nthreads == 1:
        self.nthreads=1
    else:
        if nthreads <= multiprocessing.cpu_count() and nthreads >1 and openmp.detect():
            self.nthreads = nthreads
        else:
            if nthreads > multiprocessing.cpu_count():
                raise Exception("Maximum number of threads is "+'{0:d}'.format(multiprocessing.cpu_count()))
            elif nthreads <= 1:
                raise Exception("Number of threads must be between 2 and {0:d}".format(multiprocessing.cpu_count()))
            else:
                raise Exception("OpenMP not enabled: do not set the nthreads parameter")

    #ここの部分を位置から描き直さないといけない
    #トランジットタイプは一旦省く
    self.ds_12, self.ds_13, self.ds_23 = _rsky._rsky(self.t, self.t0_1, self.per_1, self.a_1, self.inc_1, self.ecc_1, self.w_1, \
                                                     self.t0_2, self.per_2, self.a_2, self.inc_2, self.ecc_2, self.w_2, self.nthreads)


    def light_curv(self):
        #primaryかsecondaryかの違いをきちんと区別しないと厳しそう

        #3体のライトカーブを同時にだそう
        if self.limb_dark_1 == "quadratic" and self.limb_dark_2 == "quadratic":
            lc = _quadratic_ld._quadratic_ld(self.ds_12, self.ds_13, self.ds_23, self.rp_1, self.rp_2, self.u_1[0], self.u_1[1], self.u_2[0], self.u_2[1], self.nthreads)
        return lc


def get_rsky(t, t0_1, per_1, a_1, inc_1, ecc_1, w_1, t0_2, per_2, a_2, inc_2, ecc_2, w_2, nthreads):
    #https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/keplerOrbitAPI.html
    # tau_1 = ??
    # tau_2 = ??
    #軌道を計算
    #共通重心()=1体目)を持つと仮定
    ke_1 = pyasl.KeplerEllipse(a_1, per_1, e=ecc_1, tau=tau_1, omega=w_1, i=inc_1)
    ke_2 = pyasl.KeplerEllipse(a_2, per_2, e=ecc_2, tau=tau_2, omega=w_2, i=inc_2)
    #rskyのうち1体目とのものはすぐに出せる
    # ds_12 = ke_1.projPlaStDist(t)
    # ds_13 = ke_2.projPlaStDist(t)
    # ds_23 = ..
    #軌道の図を見なおしましょう。そしてxyx座標とかを自分で立式
