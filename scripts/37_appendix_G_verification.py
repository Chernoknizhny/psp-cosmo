# -*- coding: utf-8 -*-
"""
Приложение Г (v1.1): Импеданс тора и масштабная инвариантность — верификация.
Автор проверки: независимый пересчёт (рецензентский скрипт).

Проверяемые утверждения:
  A. Phi_tor = GM/(Rc^2) = 8.535 из первых принципов
  B. gamma(0.9c) = 2.294
  C. Три канала: Z_g = sqrt(Phi), Z_k = 1-beta, Z_r = 1/gamma
  D. Z = 3.457, M0 = 1/Z = 0.28924
  E. Таблица 2 (степени Phi), таблица 7 (альтернативы), таблица 10 (пределы)
  F. Декомпозиция: без вращения M0 = 0.203, Delta_M = 0.086, Delta_t = 14.1 млрд лет
  G. Порядко-инвариантное разложение через Delta_Z: 61.5% / 38.5%
  H. Frame-dragging: 2*Phi = 17.1, v_LT = 15.35c
  I. beta_coeff = e0*sqrt(Phi) = 0.3505
  J. Фазовое окно delta = 0.067 -> угловой размер 24 град.
  K. lambda_квант ~ 1e-60
  L. Масштабная инвариантность (таблица 13): Phi, Z, M0 инвариантны при R->lambda R, M->lambda M
"""
import numpy as np

G = 6.674e-8          # гравитационная постоянная, см^3/(г с^2)
c = 2.99792458e10     # скорость света, см/с
R_tor = 2.61e25       # большой радиус тора, см
M_tor = 3.0e54        # масса тора, г
T_cycle = 164.0       # длительность цикла, млрд лет
e0 = 0.12             # эксцентриситет СМЧД

results = []
def check(name, claimed, computed, tol_rel=5e-3):
    rel = abs(computed - claimed) / max(abs(claimed), 1e-300)
    status = "PASS" if rel <= tol_rel else "FAIL"
    results.append((name, claimed, computed, rel, status))

# --- A. Гравитационный потенциал тора ---
Phi = G * M_tor / (R_tor * c**2)
check("Phi_tor = GM/(Rc^2)", 8.535, Phi)

# --- B. Лоренц-фактор ---
beta = 0.9
gamma = 1 / np.sqrt(1 - beta**2)
check("gamma(0.9c)", 2.294, gamma)

# --- C. Три канала импеданса ---
Zg = np.sqrt(Phi)
Zk = 1 - beta
Zr = 1 / gamma
check("Z_g = sqrt(Phi)", 2.921, Zg)
check("Z_k = 1 - beta", 0.100, Zk)
check("Z_r = 1/gamma", 0.436, Zr)
Z = Zg + Zk + Zr
M0 = 1 / Z
check("Z = Z_g+Z_k+Z_r", 3.457, Z)
check("M0 = 1/Z", 0.28924, M0)

# --- E1. Таблица 2: проверка степени Phi ---
for p, z_cl, m_cl in [(1/3, 2.580, 0.388), (1/2, 3.457, 0.289),
                      (1, 9.066, 0.110), (2, 73.30, 0.014)]:
    Zp = Phi**p + Zk + Zr
    check(f"Табл.2 Phi^{p:.3f}: Z", z_cl, Zp)
    check(f"Табл.2 Phi^{p:.3f}: M0", m_cl, 1/Zp)

# --- E2. Таблица 7: альтернативные модели ---
Z_mult = Zg * Zk * Zr
Z_par = 1 / (1/Zg + 1/Zk + 1/Zr)
check("Мультипликативная Z", 0.127, Z_mult)
check("Мультипликативная M0", 7.86, 1/Z_mult)
check("Параллельная Z", 0.079, Z_par)
check("Параллельная M0", 12.6, 1/Z_par)

# --- E3. Таблица 10: предельные случаи ---
check("v->c: Z = sqrt(Phi)", 2.921, Zg)
check("v->c: M0", 0.342, 1 / Zg)
check("Phi=0: Z", 0.536, Zk + Zr)
check("Phi=0: M0", 1.87, 1 / (Zk + Zr))

# --- F. Декомпозиция гравитация + вращение ---
Z0 = Zg + 1 + 1                 # без вращения: Z_k = Z_r = 1
M0_v0 = 1 / Z0
t_v0 = M0_v0 * T_cycle
check("Z^(v=0)", 4.921, Z0)
check("M0^(v=0)", 0.203, M0_v0)
check("t_v0 = 33.3 млрд лет", 33.3, t_v0)
check("t_v0 vs 33.64 (время остывания)", 33.64, t_v0, tol_rel=0.015)
dM = M0 - M0_v0
dt = dM * T_cycle
check("Delta_M", 0.086, dM)
check("Delta_t = 14.1 млрд лет", 14.1, dt)
check("Delta_t vs 13.8 (2.3%)", 13.8, dt, tol_rel=0.03)
check("t_now = 47.43 млрд лет", 47.43, M0 * T_cycle)

# --- G. Порядко-инвариантная декомпозиция через Delta_Z ---
dZ = Z0 - Z
dZk, dZr = 1 - Zk, 1 - Zr
check("Delta_Z", 1.464, dZ)
check("Доля кинематического канала 61.5%", 0.615, dZk / dZ)
check("Доля релятивистского канала 38.5%", 0.385, dZr / dZ)
# Таблица 9: разложение через Delta_M (зависит от порядка) — справочно
f_k_ord1 = (1/(Z0 - dZk) - M0_v0) / (M0 - M0_v0)
f_k_ord2 = (M0 - 1/(Z0 - dZr)) / (M0 - M0_v0)
check("Табл.9 порядок k->r: 52.8%", 0.528, f_k_ord1)
check("Табл.9 порядок r->k: 69.1%", 0.691, f_k_ord2, tol_rel=0.01)

# --- H. Frame-dragging ---
check("Omega_LT/Omega_tor = 2*Phi = 17.1", 17.1, 2 * Phi)
check("v_LT = 2*Phi*beta = 15.35c", 15.35, 2 * Phi * beta)

# --- I. Коэффициент из Приложения Б ---
check("beta_coeff = e0*sqrt(Phi) = 0.3505", 0.3505, e0 * np.sqrt(Phi))

# --- J. Фазовое окно и угловой размер ---
check("delta = 0.067 -> 24 град.", 24.0, 0.067 * 360, tol_rel=0.02)

# --- K. Квантовый предел масштабирования ---
lam_q = np.sqrt(1e-121 / 0.067)
check("lambda_квант ~ 1e-60", 1e-60, lam_q, tol_rel=0.5)

# --- L. Масштабная инвариантность (таблица 13) ---
inv_ok = True
for lam in [1e-3, 1e-2, 0.1, 1, 10, 100]:
    Phi_l = G * (M_tor * lam) / ((R_tor * lam) * c**2)
    Z_l = np.sqrt(Phi_l) + (1 - beta) + 1 / gamma
    if abs(Phi_l - Phi) > 1e-6 or abs(1/Z_l - M0) > 1e-6:
        inv_ok = False
check("Табл.13: Phi, Z, M0 инвариантны по lambda", 1, int(inv_ok), tol_rel=1e-9)

# --- Отчёт ---
print("=" * 100)
print(f"{'Статус':6} | {'Проверка':42} | {'Заявлено':>12} | {'Расчёт':>14} | {'Откл.':>10}")
print("-" * 100)
for name, cl, comp, rel, st in results:
    print(f"{st:6} | {name:42} | {cl:12.6g} | {comp:14.6g} | {rel:10.2e}")
npass = sum(1 for r in results if r[4] == "PASS")
print("-" * 100)
print(f"ИТОГ: PASS {npass} / {len(results)}")
print()
print("Примечания (воспроизводимость):")
print(" 1. Формула (1) во введении напечатана с корнем над ВСЕЙ суммой: 1/sqrt(Phi+(1-b)+1/g)")
print(f"    даёт M0 = {1/np.sqrt(Phi+Zk+Zr):.4f} вместо 0.28924. Правильная запись — формулы (2)-(3).")
print(" 2. Заявленная лямбда-поправка ~0.25% при Lambda != 0 не воспроизводится:")
print(f"    (R_tor/1/sqrt(Lambda))^2 = {(R_tor/9.535e27)**2:.1e}, (10*R_tor/R_L)^2 = {(10*R_tor/9.535e27)**2:.1e}.")
print("    Требуется вывести, откуда 0.25%.")
print(" 3. Вывод Z_g = sqrt(Phi) — пропорциональность (T_cycle ∝ R_tor постулируется),")
print("    а численная валидация — согласие 1% с M_struct=0.205 и 2.3% с 13.8 млрд лет.")
print(" 4. Таблица 6: в строке lapse (ADM) значение '900' — артефакт распознавания.")
print(" 5. В п. 8.3 опечатка: 'Z_r = 1/gamma^2' вместо 'Z_r = 1/gamma'.")
