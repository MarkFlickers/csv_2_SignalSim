#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
csv_to_signalsim_traj_pwconst_v2.py

CSV (one row per epoch) -> SignalSim trajectory JSON, with near-exact 1 Hz reconstruction.

Why this variant:
Some SignalSim builds appear to interpret:
  - initVelocity.course as RADIANS (regardless of "angleUnit")
  - HorizontalTurn.angle in trajectoryList as DEGREES
This script supports that "mixed" behavior via flags.

CSV format (no header by default):
  lat, lon[, alt]
or
  lon, lat[, alt]  (use --order lonlat)

Reconstruction approach (per each dt interval):
  1) Tiny HorizontalTurn (eps_turn) to set new course
  2) Tiny ConstAcc (eps_acc) to set new speed
  3) Const for the remaining time

This makes displacement during turn/acc negligible, so 1 Hz positions match the CSV closely.

Examples:
  # Standard (everything in degrees)
  python csv_to_signalsim_traj_pwconst_v2.py track.csv out.json

  # Mixed: init course in radians, turn angles in degrees (what you observed)
  python csv_to_signalsim_traj_pwconst_v2.py track.csv out.json --init-course-unit rad --turn-angle-unit degree --angleunit-field degree

  # Use a template JSON (only replace "trajectory")
  python csv_to_signalsim_traj_pwconst_v2.py track.csv out.json --template base.json --init-course-unit rad
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

a_wgs84 = 6378137.0
e2_wgs84 = 0.00669437999014

def R_ecef_to_enu(lat0: float, lon0: float) -> np.ndarray:
    """ Матрица поворота ECEF->ENU в точке (lat0, lon0). """
    sin_lat, cos_lat = np.sin(lat0), np.cos(lat0)
    sin_lon, cos_lon = np.sin(lon0), np.cos(lon0)

    return np.array([
        [-sin_lon,            cos_lon,           0.0],
        [-sin_lat*cos_lon, -sin_lat*sin_lon,  cos_lat],
        [ cos_lat*cos_lon,  cos_lat*sin_lon,  sin_lat]
    ])

def llh_to_ecef(lat: float, lon: float, alt: float) -> Tuple[float, float, float]:
    """LLH (рад, рад, м) -> ECEF (м)."""
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    sin_lon, cos_lon = np.sin(lon), np.cos(lon)

    N = a_wgs84 / np.sqrt(1.0 - e2_wgs84 * sin_lat**2)

    x = (N + alt) * cos_lat * cos_lon
    y = (N + alt) * cos_lat * sin_lon
    z = (N * (1.0 - e2_wgs84) + alt) * sin_lat
    return x, y, z

def ecef_delta_to_enu(dx: float, dy: float, dz: float, lat0: float, lon0: float) -> Tuple[float, float, float]:
    """ECEF-вектор (приращение) -> ENU-вектор в точке (lat0, lon0)."""
    R = R_ecef_to_enu(lat0, lon0)
    e, n, u = R @ np.array([dx, dy, dz])
    return float(e), float(n), float(u)

def wrap180(deg: float) -> float:
    return ((deg + 180.0) % 360.0) - 180.0

prev_ang = 0.0
def course_deg(east_m: float, north_m: float) -> float:
    """Course clockwise from north, degrees [0, 360)."""
    global prev_ang
    if((east_m == 0) and (north_m == 0)):
        ang = prev_ang
    else:
        ang = math.degrees(math.atan2(east_m, north_m))
    prev_ang = ang
    return (ang + 360.0) % 360.0


def read_track_csv(path: Path, order: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(path, header=None)
    if df.shape[1] < 2:
        raise ValueError("CSV must have at least 2 columns (lat, lon) or (lon, lat).")

    c0 = df.iloc[:, 0].astype(float).to_numpy()
    c1 = df.iloc[:, 1].astype(float).to_numpy()
    alt = df.iloc[:, 2].astype(float).to_numpy() if df.shape[1] >= 3 else None

    if order == "latlon":
        lat, lon = c0, c1
    elif order == "lonlat":
        lon, lat = c0, c1
    else:
        raise ValueError("order must be 'latlon' or 'lonlat'.")
    return lat, lon, alt


def estimate_speed_course_vert_speed(
    lat_deg: np.ndarray,
    lon_deg: np.ndarray,
    alt: np.ndarray,
    dt: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    if alt is None:
        alt = np.zeros_like(lat)
        no_vertical = True
    else:
        no_vertical = False

    speed_hor = []
    course = []
    speed_vert = []

    for i in range(len(lat) - 1):
        x1, y1, z1 = llh_to_ecef(lat[i], lon[i], alt[i])
        x2, y2, z2 = llh_to_ecef(lat[i + 1], lon[i + 1], alt[i + 1])

        dx, dy, dz = x2 - x1, y2 - y1, z2 - z1

        lat_ref = 0.5 * (lat[i] + lat[i + 1])
        lon_ref = 0.5 * (lon[i] + lon[i + 1])

        e, n, u = ecef_delta_to_enu(dx, dy, dz, lat_ref, lon_ref)

        dist_hor = np.hypot(e, n)
        speed_hor.append(dist_hor / dt)
        course.append(course_deg(e, n))
        speed_vert.append(u / dt)

    if no_vertical == True:
        speed_vert = np.zeros_like(speed_vert)
    return np.asarray(speed_hor), np.asarray(course), np.asarray(speed_vert)


def build_pwconst_segments(hor_speed: np.ndarray,
                           course: np.ndarray,
                           vert_speed: np.ndarray,
                           dt: float,
                           eps_turn: float,
                           eps_acc: float) -> List[Dict[str, Any]]:
    if len(hor_speed) == 0:
        return [{"type": "Const", "time": 0.0}]

    segs: List[Dict[str, Any]] = []
    cur_hor_speed = float(hor_speed[0])
    cur_course = float(course[0])
    cur_vert_speed = float(vert_speed[0])
    segs.append({"type": "Const", "time": float(dt)})

    for i in range(1, len(hor_speed)):
        tgt_hor_speed = float(hor_speed[i])
        tgt_course = float(course[i])
        tgt_vert_speed = float(vert_speed[i])
        cur_course_speed = np.hypot(cur_hor_speed, cur_vert_speed)
        if cur_hor_speed > 1e-9:
            tgt_course_speed = cur_course_speed * (tgt_hor_speed/cur_hor_speed)
        else:
            tgt_course_speed = cur_course_speed

        used = 0.0

        d_course = wrap180(tgt_course - cur_course)

        if abs(d_course) > 1e-8 and eps_turn > 0:
            segs.append({"type": "HorizontalTurn", "time": float(eps_turn), "angle": float(d_course)})
            cur_course = (cur_course + d_course) % 360.0
            used += eps_turn

        if abs(tgt_course_speed - cur_course_speed) > 1e-12 and eps_acc > 0:
            segs.append({"type": "ConstAcc", "time": float(eps_acc), "speed": float(tgt_course_speed)})
            if cur_course_speed > 1e-12:
                k = tgt_course_speed / cur_course_speed
                cur_hor_speed *= k
                cur_vert_speed *= k
            else:
                cur_hor_speed = tgt_hor_speed
            used += eps_acc

        if abs(tgt_vert_speed - cur_vert_speed) > 1e-12 and eps_acc > 0:
            segs.append({"type": "VerticalAcc", "time": float(eps_acc), "speed": float(tgt_vert_speed)})
            cur_vert_speed = tgt_vert_speed
            used += eps_acc

        t_rem = dt - used
        if t_rem < 0:
            t_rem = 0.0
        if t_rem > 0:
            segs.append({"type": "Const", "time": float(t_rem)})

    return segs


def convert_course(value_deg: float, unit: str) -> float:
    if unit == "degree":
        return float(value_deg)
    if unit == "rad":
        return float(value_deg) * (math.pi / 180.0)
    raise ValueError("init-course-unit must be 'degree' or 'rad'")


def convert_turn_angles(segs: List[Dict[str, Any]], unit: str) -> None:
    if unit == "degree":
        return
    if unit == "rad":
        k = math.pi / 180.0
        for seg in segs:
            if seg.get("type") == "HorizontalTurn":
                seg["angle"] = float(seg["angle"]) * k
        return
    raise ValueError("turn-angle-unit must be 'degree' or 'rad'")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path, help="Input CSV track (1 row per epoch).")
    ap.add_argument("json_out", type=Path, help="Output JSON.")
    ap.add_argument("--order", choices=["latlon", "lonlat"], default="latlon")
    ap.add_argument("--dt", type=float, default=1.0, help="Epoch spacing in seconds.")
    ap.add_argument("--eps-turn", type=float, default=0.001, help="Tiny HorizontalTurn duration (s).")
    ap.add_argument("--eps-acc", type=float, default=0.001, help="Tiny ConstAcc duration (s).")
    ap.add_argument("--init-course-unit", choices=["degree", "rad"], default="degree")
    ap.add_argument("--turn-angle-unit", choices=["degree", "rad"], default="degree")
    ap.add_argument("--angleunit-field", choices=["degree", "rad"], default="degree")
    ap.add_argument("--name", type=str, default="csv_trajectory_pwconst")
    ap.add_argument("--template", type=Path, default=None)

    args = ap.parse_args(argv)

    if args.eps_turn + 2 * args.eps_acc >= args.dt:
        raise ValueError("eps_turn + 2*eps_acc must be < dt.")

    lat, lon, alt = read_track_csv(args.csv, args.order)
    speed, course, vertical_speed = estimate_speed_course_vert_speed(lat, lon, alt, args.dt)
    segs = build_pwconst_segments(speed, course, vertical_speed, args.dt, args.eps_turn, args.eps_acc)

    # Apply requested unit for turn angles
    convert_turn_angles(segs, args.turn_angle_unit)

    # Init course: compute from first interval in degrees, then convert if requested
    init_course_deg = float(course[0]) if len(course) else 0.0
    init_course = convert_course(init_course_deg, args.init_course_unit)
    init_vert_speed = vertical_speed[0]

    trajectory = {
        "name": args.name,
        "initPosition": {
            "type": "LLA",
            "format": "d",
            "longitude": float(lon[0]),
            "latitude": float(lat[0]),
            "altitude": float(0 if alt is None else alt[0]),
        },
        "initVelocity": {
            "type": "SCU",
            "speedUnit": "mps",
            "angleUnit": args.angleunit_field,
            "speed": float(speed[0]) if len(speed) else 0.0,
            "course": float(init_course),
            "up": float(init_vert_speed),
        },
        "trajectoryList": segs,
    }

    if args.template:
        cfg = json.loads(args.template.read_text(encoding="utf-8"))
        cfg["trajectory"] = trajectory
    else:
        cfg = {
            "version": 1.0,
            "description": "trajectory generated from CSV (pw-const reconstruction)",
            "time": {"type": "UTC", "year": 2026, "month": 1, "day": 1, "hour": 0, "minute": 0, "second": 0},
            "trajectory": trajectory,
        }

    args.json_out.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote: {args.json_out}")


if __name__ == "__main__":
    main()
