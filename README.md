# csv_to_signalsim_traj.py — генерация траектории SignalSim из CSV

Скрипт конвертирует CSV-трек (по одной строке на эпоху) в JSON-описание траектории для SignalSim (формат `trajectory`).
Основная цель — **максимально точно воспроизвести исходные точки на частоте 1 Гц** при использовании сегментной модели SignalSim.

Ключевая особенность: поддерживается «смешанный» режим единиц угла, который встречается в некоторых сборках SignalSim:
- `initVelocity.course` ожидается **в радианах** (даже если `angleUnit` в JSON стоит `degree`);
- `HorizontalTurn.angle` в `trajectoryList` ожидается **в градусах**.

Скрипт позволяет явно выбрать единицы для обоих полей.

---

## Зависимости

- Python 3.9+
- `numpy`
- `pandas`

Установка:
```bash
pip install numpy pandas
```

---

## Входной CSV

По умолчанию — без заголовка, одна строка на эпоху:

- `lat, lon[, alt]`  (широта, долгота, высота)
- либо `lon, lat[, alt]` (если задано `--order lonlat`)

Единицы:
- lat/lon — **в градусах**
- alt — **в метрах** (если нет столбца высоты, берётся 0)

Интервал между строками задаётся параметром `--dt` (по умолчанию 1.0 с).

---

## Как устроена генерация траектории (pw-const)

SignalSim работает сегментами. Чтобы 1-Гц точки совпадали с CSV, на каждый интервал `dt` формируется последовательность:

1. очень короткий поворот `HorizontalTurn` длительностью `eps_turn`  
2. очень короткое изменение скорости `ConstAcc` длительностью `eps_acc`  
3. движение по прямой `Const` на оставшееся время

За счёт малых `eps_*` смещение во время поворота/разгона пренебрежимо мало, и траектория почти проходит через исходные точки.

По умолчанию:
- `eps_turn = 0.001` с
- `eps_acc  = 0.001` с

Важно: должно выполняться `eps_turn + eps_acc < dt`.

---

## Базовый запуск

Создать JSON, содержащий только `time` и `trajectory`:
```bash
python csv_to_signalsim_traj.py track.csv out.json
```

---

## Работа через шаблон (рекомендуется)

Если у вас уже есть полноценный JSON-конфиг SignalSim (ephemeris/output/power/маски и т.п.), удобнее использовать его как шаблон.
Скрипт заменит **только** поле `trajectory`.

```bash
python csv_to_signalsim_traj.py track.csv Out.json --template base.json
```

---

## Рабочий режим (init course в радианах, повороты в градусах)

```bash
python csv_to_signalsim_traj.py track.csv Out.json ^
  --template out.json ^
  --init-course-unit rad ^
  --turn-angle-unit degree ^
  --angleunit-field degree
```

Пояснения:
- `--init-course-unit rad` — записать `initVelocity.course` в радианах;
- `--turn-angle-unit degree` — оставить `HorizontalTurn.angle` в градусах;
- `--angleunit-field degree` — что писать в `initVelocity.angleUnit` (некоторые сборки это поле игнорируют, но для совместимости можно оставить `degree`).

---

## Параметры командной строки

- `--order {latlon,lonlat}` — порядок столбцов в CSV (по умолчанию `latlon`)
- `--dt <сек>` — шаг по времени между строками CSV (по умолчанию `1.0`)
- `--eps-turn <сек>` — длительность «микроповорота» (по умолчанию `0.001`)
- `--eps-acc <сек>` — длительность «микроразгона» (по умолчанию `0.001`)
- `--init-course-unit {degree,rad}` — единицы для `initVelocity.course`
- `--turn-angle-unit {degree,rad}` — единицы для `HorizontalTurn.angle`
- `--angleunit-field {degree,rad}` — значение поля `initVelocity.angleUnit`, записываемое в JSON
- `--template <json>` — шаблон JSON, в котором будет заменён только раздел `trajectory`
- `--name <строка>` — имя траектории `trajectory.name` (по умолчанию `csv_trajectory_pwconst`)


