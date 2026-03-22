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

EARTH_R = 6371000.0  # meters


def wrap180(deg: float) -> float:
    return ((deg + 180.0) % 360.0) - 180.0


def course_deg(east_m: float, north_m: float) -> float:
    """Course clockwise from north, degrees [0, 360)."""
    ang = math.degrees(math.atan2(east_m, north_m))
    return (ang + 360.0) % 360.0


def read_track_csv(path: Path, order: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(path, header=None)
    if df.shape[1] < 2:
        raise ValueError("CSV must have at least 2 columns (lat, lon) or (lon, lat).")

    c0 = df.iloc[:, 0].astype(float).to_numpy()
    c1 = df.iloc[:, 1].astype(float).to_numpy()
    alt = df.iloc[:, 2].astype(float).to_numpy() if df.shape[1] >= 3 else np.zeros_like(c0)

    if order == "latlon":
        lat, lon = c0, c1
    elif order == "lonlat":
        lon, lat = c0, c1
    else:
        raise ValueError("order must be 'latlon' or 'lonlat'.")
    return lat, lon, alt


def estimate_speed_course(lat_deg: np.ndarray, lon_deg: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Approximate EN displacements on a tangent plane:
      north ≈ dlat * R
      east  ≈ dlon * R * cos(lat_mid)
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)

    dlat = np.diff(lat)
    dlon = np.diff(lon)
    latm = 0.5 * (lat[:-1] + lat[1:])

    north = dlat * EARTH_R
    east = dlon * EARTH_R * np.cos(latm)

    dist = np.sqrt(north * north + east * east)
    speed = dist / dt
    course = np.array([course_deg(e, n) for e, n in zip(east, north)], dtype=float)
    return speed, course


def build_pwconst_segments(speed: np.ndarray,
                           course: np.ndarray,
                           dt: float,
                           eps_turn: float,
                           eps_acc: float) -> List[Dict[str, Any]]:
    if len(speed) == 0:
        return [{"type": "Const", "time": 0.0}]

    segs: List[Dict[str, Any]] = []
    cur_speed = float(speed[0])
    cur_course = float(course[0])

    for i in range(len(speed)):
        tgt_speed = float(speed[i])
        tgt_course = float(course[i])

        used = 0.0

        d_course = wrap180(tgt_course - cur_course)
        if abs(d_course) > 1e-12 and eps_turn > 0:
            segs.append({"type": "HorizontalTurn", "time": float(eps_turn), "angle": float(d_course)})
            cur_course = (cur_course + d_course) % 360.0
            used += eps_turn

        if abs(tgt_speed - cur_speed) > 1e-12 and eps_acc > 0:
            segs.append({"type": "ConstAcc", "time": float(eps_acc), "speed": float(tgt_speed)})
            cur_speed = tgt_speed
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

    if args.eps_turn + args.eps_acc >= args.dt:
        raise ValueError("eps_turn + eps_acc must be < dt. Reduce eps values.")

    lat, lon, alt = read_track_csv(args.csv, args.order)
    speed, course = estimate_speed_course(lat, lon, args.dt)
    segs = build_pwconst_segments(speed, course, args.dt, args.eps_turn, args.eps_acc)

    # Apply requested unit for turn angles
    convert_turn_angles(segs, args.turn_angle_unit)

    # Init course: compute from first interval in degrees, then convert if requested
    init_course_deg = float(course[0]) if len(course) else 0.0
    init_course = convert_course(init_course_deg, args.init_course_unit)

    trajectory = {
        "name": args.name,
        "initPosition": {
            "type": "LLA",
            "format": "d",
            "longitude": float(lon[0]),
            "latitude": float(lat[0]),
            "altitude": float(alt[0]),
        },
        "initVelocity": {
            "type": "SCU",
            "speedUnit": "mps",
            "angleUnit": args.angleunit_field,
            "speed": float(speed[0]) if len(speed) else 0.0,
            "course": float(init_course),
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
