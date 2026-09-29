# psp-cosmo

Репозиторий с кодом для сравнения моделей PSP и ΛCDM.

## Как запустить

```bash
git clone https://github.com/Chernoknizhny/psp-cosmo.git
cd psp-cosmo
python -m venv venv
# Windows: venv\Scripts\activate
# Linux/macOS: source venv/bin/activate
pip install -r requirements.txt
dvc pull
Скрипты
python scripts/compute_chi2_desi.py — посчитать χ² на DESI DR1.
python scripts/plot_alpha_H_scan.py — построить график Δχ² vs α_H.
python scripts/retro_forecast.py — полный расчёт.
Структура
src/psp_model.py — формула H(z) для PSP.
src/lcdm_model.py — формула H(z) для ΛCDM.
src/stats.py — функции для χ², BIC, Δχ².
config/params.yaml — параметры B=2.591, w=0.0367 (фиксированы).
scripts/ — скрипты запуска.
data/ — данные (через DVC).
results/ — результаты.
Результаты (числа)
Модель	χ² (DESI DR1)	Примечание
PSP	3.77	B, w не фитились
ΛCDM	8.56	w = -1
Δχ²	-4.79	PSP лучше
Данные
Pantheon+: Scolnic et al. (2022).
DESI DR1: Aghamousa et al. (2022).
CMB: Planck Collaboration.
Лицензия
MIT.
