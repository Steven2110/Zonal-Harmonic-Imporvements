# Files path
SUBFOLDERS = ['x_1', 'x_2', 'x_3', 'vx_1', 'vx_2', 'vx_3', 'J20']
PATH_VARIATION_LOG = "log/variation_selection.log"
PATH_ORBIT_LOG = "log/orbit_improvements.log"
PATH_OUT = "EPH.OUT"
PATH_IN = "iszm_puc.in"
PATH_EXE = "iszm_puc.exe"
PATH_360 = "MODTOOLS/garm360.in"
PATH_360_SOURCE = "MODTOOLS/originalgarm360.in"
PATH_FIRST_STEP = "data/variation_selection/first_step_data.csv"
PATH_VARIATION_STEP = "data/variation_selection/{}/variation_step_data_{}.csv"
PATH_OBSERVATION_DATA = "data/orbit_improvements/observation_data.csv"
PATH_CALCULATION_DATA = "data/orbit_improvements/calculation_data_iteration_{}.csv"
PATH_RESULT = "data/orbit_improvements/result.csv"
PATH_FIGURE = "pics/derivative_variation_{}.png"

# J20 Coefficient
J20 = -0.186987635955E-09

# Arcsecond
ARCSECOND = 0.01

# Date
T0 = 2458369.5

# Gravitational parameters for Earth
GM = 398600.4418

# Orbital elements for Geostationary Satellite
SEMI_MAJOR_AXIS = 42165
ECCENTRICITY = 0.001
INCLINATION = 1
LONGITUDE_OF_ASCENDING_NODE = 0
ARGUMENT_OF_PERIAPSIS = 0
MEAN_ANOMALY = 0
T = 86400
N = 100

INPUT_TEMPLATE = """1   РЕЖИМ (1 - прогноз; 2 - улучшение орбиты)
2018  9  8  0  0  0.000   Начальная эпоха (TT)
 1     Число спутников
  {}  {}  {}
  {}  {}  {}
  
         ПРОГНОЗ
2018  9  8  0  0  0.000  Начальный момент прогноза (TT)
{}  {}  {}  {}  {}  {}  Конечный момент прогноза (TT) 
 {} Шаг выдачи (сек)
1E-3   Ошибка большой полуоси (км)

УЛУЧШЕНИЕ
obsvyb.in  Файл с наблюдениями
5           Число наблюдений
2822.172315  2201.434565  5279.171158  координаты обсерватории (км)
0.10         Точность улучшения (сходимости) (км)
1E-5 1E-8    Вариации начальных координат и скоростей (км и км/c)
0.1        Множитель, улучшающий сходимость
2            Начальные условия (1 - из файла; 2 - вычисляются в программе)
8000        Приблизительный радиус-вектор (км) на момент 1-го наблюдения

ИНТЕГРИРОВАНИЕ
   10.  Постоянный шаг интегрирования (сек; для отриц. параметра)
   19   Порядок интегратора (от 7 до 39 через 4)
   10   Параметр интегратора
 1000   Интервал промежуточных выдач на экран (в шагах интегрирования)

ВОЗМУЩЕНИЯ
 2  0   Гармоники геопотенциала (NM)
    0   Луна
    0   Солнце
    0   Световое давление и ПР эффект
0 0 0   Релятивистские эффекты (Ф_0, Ф_1, Ф_2)
    0   Приливы
    0   Атмосфера
  100.  Высота сгорания (км)

СПУТНИК
  500.   Масса (кг)
  0.5   Площадь миделева сечения (м^2)
   2.   Коэф-т лобового сопротивления
   2.   Коэф-т отражения
"""

GARM_TEMPLATE = "   2   0 {}  0.000000000000E+00\n"