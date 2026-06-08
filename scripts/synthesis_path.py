"""
RxnFlow 학습 로그(SQLite)에서 합성 경로(traj)를 꺼내,
dual_target_hits.csv 의 각 hit 에 합성 경로 초안을 붙인다.

traj 는 JSON 문자열: [{"step":0,"smiles":...,"action":["FirstBlock","<block>"]},
                      {"step":1,...,"action":["BiRxn","<template>","<block>"]}, ...]

사용:
  python ./scripts/synthesis_path.py \
    --train_dir ./logs/pdl1_6vqn/train \
    --hits_csv ./logs/dual_target_hits.csv \
    --out ./logs/dual_target_hits_with_path.csv
"""

import argparse
import glob
import json
import sqlite3
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def canonical(smi: str) -> str | None:
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None


def load_traj_map(train_dir: str) -> dict[str, str]:
    """SQLite 로그들에서 canonical_smiles -> traj(JSON) 매핑 구축."""
    dbs = sorted(glob.glob(str(Path(train_dir) / "generated_objs_*.db")))
    if not dbs:
        raise FileNotFoundError(f"generated_objs_*.db 없음: {train_dir}")
    smi2traj: dict[str, str] = {}
    for db in dbs:
        con = sqlite3.connect(db)
        try:
            # 테이블명은 보통 'results' (gflownet SQLiteLog 기본)
            tbls = [r[0] for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
            tbl = "results" if "results" in tbls else tbls[0]
            for smi, traj in con.execute(f"SELECT smi, traj FROM {tbl}"):
                if not smi or not traj:
                    continue
                c = canonical(smi)
                if c and c not in smi2traj:
                    smi2traj[c] = traj
        finally:
            con.close()
    print(f"로그에서 합성 경로 {len(smi2traj)}개 로드")
    return smi2traj


def summarize(traj_json: str) -> str:
    """traj JSON -> 사람이 읽는 한 줄 요약 (블록/반응 순서)."""
    try:
        steps = json.loads(traj_json)
    except Exception:
        return ""
    parts = []
    for s in steps:
        a = s.get("action", [])
        if not a:
            continue
        if a[0] == "FirstBlock":
            parts.append(f"start:{a[1]}")
        elif a[0] == "UniRxn":
            parts.append(f"uni[{a[1]}]")
        elif a[0] == "BiRxn":
            parts.append(f"+{a[2]} via[{a[1]}]")
        elif a[0] == "Stop":
            parts.append("stop")
    return "  ->  ".join(parts)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--train_dir", required=True, help="예: ./logs/pdl1_6vqn/train")
    p.add_argument("--hits_csv", required=True, help="dual_target_screen.py 결과 csv")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    smi2traj = load_traj_map(args.train_dir)
    df = pd.read_csv(args.hits_csv)

    smi_col = "smiles" if "smiles" in df.columns else df.columns[1]
    df["_canon"] = df[smi_col].map(canonical)
    df["synthesis_path_json"] = df["_canon"].map(lambda c: smi2traj.get(c, ""))
    df["synthesis_path"] = df["synthesis_path_json"].map(summarize)
    matched = (df["synthesis_path_json"] != "").sum()
    df = df.drop(columns=["_canon"])

    df.to_csv(args.out, index=False)
    print(f"저장: {args.out}")
    print(f"  hit {len(df)}개 중 합성 경로 매칭: {matched}개")
    if matched < len(df):
        print("  (매칭 안 된 건 PD-L1 생성 로그에 그 SMILES의 traj가 없거나 canonical 불일치)")


if __name__ == "__main__":
    main()