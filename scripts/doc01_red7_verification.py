#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
doc01_red7_verification.py — независимая проверка чисел ДОКУМЕНТА_01_PSP_ред7
(сводный обзор серии PSP, «честный протокол»).

Все константы зашиты из документа №01 (ред. 7) и ссылочных документов
(№03, №11, №16, №35, №36, №41). Свободных подгоночных параметров нет.
Запуск:  python doc01_red7_verification.py
Зависимости: numpy, scipy
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# ----------------------------------------------------------------------
# Константы модели (источник: №01 ред.7, разделы 2 и 4)
# ----------------------------------------------------------------------
PHI_TOR    = 8.535        # гравитационный потенциал тора [№03]
GAMMA      = 2.294        # при v = 0.9c [№03 §4.5]
DELTAM_MAX = 13.8/164.0   # фаза наблюдателя: t0/T [№16]
T_CYCLE    = 164.0        # длина цикла, млрд лет (L = 2*pi*a)
H0_PLANCK  = 67.36        # нормировочное значение [Planck 2018]
OM_PLANCK  = 0.315
H0_MODEL   = 70.9         # модельное предсказание из геометрии [№16]
OM_MODEL   = 0.30
E_C        = 0.316        # эксцентриситет оболочки: sqrt(r/R) [№16]
ALPHA_H    = 0.1464       # заявлено [№36]
BETA       = 0.3506       # заявлено [№36]
U_MAX      = 20.0         # возраст интегрируется до z = e^20 ~ 5e8
N_GRID     = 100000


def g_factor(dm):
    """g(DeltaM) = exp(-DeltaM/w) * (1 + B*DeltaM/(DeltaM + w))."""
    dm = np.asarray(dm, dtype=float)
    return np.exp(-dm/w_param)*(1.0 + B*dm/(dm + w_param))


class TimeScale:
    """Lookback-время t_lb(z) в ACDM-шкале (H0, Om); DeltaM = t_lb / T_cycle."""

    def __init__(self, H0, Om):
        self.H0, self.Om = H0, Om

    def t_lb_Gyr(self, z):
        E = lambda zz: np.sqrt(self.Om*(1.0+zz)**3 + (1.0-self.Om))
        val, _ = quad(lambda u: 1.0/E(np.expm1(u)), 0.0, np.log1p(z), limit=200)
        return val/self.H0*977.8

    def dm_of_z(self, z):
        return self.t_lb_Gyr(z)/T_CYCLE

    def w_de(self, z):
        """w_DE: rho_DE = rho_L*g  =>  w_DE = (1/3) dln g / dln(1+z)."""
        z = float(np.clip(z, 0.0, None))
        h = 1e-5
        z1 = max(z, 0.0); z2 = z + h
        dlng = (np.log(g_factor(self.dm_of_z(z2))) - np.log(g_factor(self.dm_of_z(z1))))/h
        return (1.0+z)/3.0*dlng

    def age(self, with_g=True):
        """Возраст, млрд лет: интегрирование по u = ln(1+z) на плотной сетке.
        t_lb на сетке считается накопительной суммой (не через quad) — это
        устойчиво и быстро; точность трапеций на 1e5 точках >> 1e-6."""
        u = np.linspace(0.0, U_MAX, N_GRID)
        z = np.expm1(u)
        Om = self.Om
        E_acdm = np.sqrt(Om*(1.0+z)**3 + (1.0-Om))
        inv = 1.0/E_acdm
        dm = (np.concatenate([[0.0], np.cumsum((inv[1:]+inv[:-1])/2.0*np.diff(u))])
              / self.H0*977.8)/T_CYCLE
        if with_g:
            E = np.sqrt(Om*(1.0+z)**3 + (1.0-Om)*g_factor(dm))
        else:
            E = E_acdm
        return np.trapezoid(1.0/E, u)*977.8/self.H0


results = []

def check(name, got, expected, tol):
    results.append((name, float(got), float(expected), float(tol), abs(got-expected) <= tol))


print("="*78)
print("ПРОВЕРКА ДОКУМЕНТА_01_PSP_ред7 (сводный обзор PSP, честный протокол)")
print("="*78)

# --- 1. Геометрические константы [ГЕО] --------------------------------
print("\n[1] Геометрические константы [ГЕО]")
B = PHI_TOR/(GAMMA+1.0)
w_param = DELTAM_MAX/GAMMA
print(f"    B = Phi_tor/(gamma+1)   = {B:.4f}   (заявлено 2.591)")
print(f"    w = DeltaM_max/gamma    = {w_param:.4f}   (заявлено 0.0367)")
print(f"    DeltaM_max = 13.8/164   = {DELTAM_MAX:.4f}   (заявлено 0.0841)")
check("B  = Phi_tor/(gamma+1)",  B, 2.591, 5e-4)
check("w  = DeltaM_max/gamma",   w_param, 0.0367, 5e-4)
check("DeltaM_max = 13.8/164",   DELTAM_MAX, 0.0841, 5e-4)

# --- 2. Свойства g(DeltaM) --------------------------------------------
print("\n[2] Свойства фактора g(DeltaM)")
xs = np.linspace(1e-6, 0.5, 200000)
vals = g_factor(xs)
i_max = int(np.argmax(vals))
dm_cross = brentq(lambda x: float(g_factor(x))-1.0, 1e-6, 0.5)
ts = TimeScale(H0_PLANCK, OM_PLANCK)
zz = np.expm1(np.linspace(0.0, U_MAX, N_GRID))
E_g = np.sqrt(OM_PLANCK*(1.0+zz)**3 + (1.0-OM_PLANCK))
inv_g = 1.0/E_g
uu = np.linspace(0.0, U_MAX, N_GRID)
tlb_grid = np.concatenate([[0.0], np.cumsum((inv_g[1:]+inv_g[:-1])/2.0*np.diff(uu))])/H0_PLANCK*977.8
dm_curve = tlb_grid/T_CYCLE
z_cross = float(np.interp(dm_cross, dm_curve, zz))
print(f"    g_max = {vals[i_max]:.4f} при DeltaM = {xs[i_max]:.5f} (заявлено 1.184 при 0.0104)")
print(f"    пересечение g=1 при DeltaM = {dm_cross:.5f} -> z = {z_cross:.2f} (заявлено z~0.42)")
check("g_max",           vals[i_max], 1.184, 5e-3)
check("DeltaM(g_max)",   xs[i_max],  0.0104, 5e-4)
check("z(g=1)",          z_cross,    0.42,   0.03)

# --- 3. Форма w_DE(z) [СЛЕП., фиксация 30.09.2026] --------------------
print("\n[3] Форма w_DE(z): rho_DE = rho_L*g(DeltaM)  [СЛЕП., фикс. 30.09.2026]")
w0 = -1.0 + ts.w_de(0.0)
z_sw = brentq(lambda z: ts.w_de(z), 0.01, 1.0)
zzs = np.linspace(0.01, 3.0, 600)
ws = -1.0 + np.array([ts.w_de(z) for z in zzs])
i_min = int(np.argmin(ws))
print(f"    w_DE(0)        = {w0:+.3f}  (заявлено +0.28)")
print(f"    переход w=-1   при z = {z_sw:.3f}  (заявлено ~0.13)")
print(f"    минимум w_DE   = {ws[i_min]:+.3f} при z = {zzs[i_min]:.2f}  (заявлено -1.41 при ~0.6)")
print(f"    w_DE(z=2)      = {-1.0+ts.w_de(2.0):+.3f}  -> «возврат к -1 при z~2» асимптотический")
check("w_DE(0)",     w0,        0.28,  0.03)
check("z(w_DE=-1)",  z_sw,      0.13,  0.02)
check("min w_DE",    ws[i_min], -1.41, 0.03)
check("z(min w_DE)", zzs[i_min], 0.6,   0.05)

# --- 4. Кросс-чек alpha_H / beta [№36] --------------------------------
print("\n[4] Кросс-чек: beta/alpha_H [№36]")
ratio = BETA/ALPHA_H
ratio_geom = PHI_TOR**0.5*(1.0-E_C**2)/(1.0+E_C**2)   # e0 сокращается
print(f"    заявленное отношение: {ratio:.3f}; из геометрии: {ratio_geom:.3f} (e={E_C})")
check("beta/alpha_H (геометрия)", ratio_geom, ratio, 0.02)

# --- 5. Возраст --------------------------------------------------------
print("\n[5] Возраст Вселенной")
age_acdm_p = ts.age(with_g=False)
age_psp_p  = ts.age(with_g=True)
ts_m = TimeScale(H0_MODEL, OM_MODEL)
age_psp_m = ts_m.age(with_g=True)
print(f"    ACDM (планковская нормировка): {age_acdm_p:.2f} млрд лет  (заявлено 13.79)")
print(f"    PSP с g (планковская):         {age_psp_p:.2f} млрд лет  (заявлено 13.83)")
print(f"    PSP с g (H0=70.9, Om=0.3):     {age_psp_m:.2f} млрд лет  (заявлено 13.1)")
check("возраст ACDM (планковская)", age_acdm_p, 13.79, 0.05)

# --- 6. Самосогласованность DeltaM(z) [открытый пункт 7.1] -------------
print("\n[6] DeltaM(2.33) и открытый пункт 7.1")
dm23_p = ts.dm_of_z(2.33)
dm23_m = ts_m.dm_of_z(2.33)
print(f"    DeltaM(2.33), ACDM-шкала (67.36/0.315): {dm23_p:.4f}")
print(f"    DeltaM(2.33), ACDM-шкала (70.9/0.300):  {dm23_m:.4f}")
print(f"    заявлено 0.1135 > DeltaM_max = 0.0841 — НЕ воспроизводится на стандартных")
print(f"    lookback-шкалах; нужно точное определение DeltaM(z) из №11 ред.2")

# ----------------------------------------------------------------------
print("\n" + "="*78)
print("СВОДКА")
print("="*78)
print(f"{'Проверка':36s} {'Получено':>9s} {'Заявлено':>9s} {'Допуск':>7s}  Статус")
print("-"*78)
for name, got, exp, tol, ok in results:
    print(f"{name:36s} {got:9.4f} {exp:9.4f} {tol:7.4f}  {'PASS' if ok else 'FAIL'}")
print("-"*78)
n_fail = sum(1 for r in results if not r[4])
print(f"Пройдено: {len(results)-n_fail}/{len(results)}")
print("""
НЕ воспроизводится в рамках документа №01 ([НЕ ПОДТВ.], см. разделы 5-6):
  * PSP-возраст 13.83 (планковская нормировка) и 13.1 (H0=70.9);
  * DeltaM(2.33) = 0.1135 > DeltaM_max (открытый пункт 7.1).
  * «Возврат w_DE к -1 при z~2» асимптотический: при z=2 отклонение ~0.2.
Для закрытия нужна таблица DeltaM(z) из №11 ред.2 (или №41).
""")
