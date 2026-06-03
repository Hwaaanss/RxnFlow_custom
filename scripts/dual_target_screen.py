"""
이중 타겟 스크리닝 (통합 스크립트, verbose 진행률 포함)
=========================================================
1) 1차 타겟 run 의 oracle SDF (docking/oracle*.sdf) 에서 SMILES + 1차 Vina 점수 추출
2) 그 후보들을 2차 타겟(임의 PDB)에 재도킹  -- 단백질/리간드/좌표를 CLI 인자로 받음
3) 두 점수를 합쳐 dual-target hit 으로 정렬 -> CSV 저장

2차 타겟은 하드코딩 없음. 네가 가진 어떤 PDB 든 --target2_protein 으로 넣으면 됨.
--verbose 를 주면 재도킹을 배치로 쪼개 진행률 / 경과시간 / ETA 를 출력.

사용 예:
  python dual_target_screen.py \
    --docking_dir ./logs/pdl1_6vqn/docking \
    --target1_name PD-L1 --target2_name TGF-bR1 \
    --target2_protein ./targets/1vjy.pdb \
    --target2_ref_ligand ./targets/ligand_1vjy.pdb \
    --out ./logs/dual_target_hits.csv \
    --pre_topk 2000 --diverse --final_topk 1000 \
    --num_workers 6 --verbose --batch 500
"""

import argparse
import glob
import time
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs, RDLogger

from rxnflow.tasks.utils.unidock import VinaReward

RDLogger.DisableLog("rdApp.*")


# ---------- 1) 1차 타겟 oracle SDF 에서 후보 추출 ----------
def extract_candidates(docking_dir: str) -> dict[str, float]:
    sdfs = sorted(glob.glob(str(Path(docking_dir) / "oracle*.sdf")))
    if not sdfs:
        raise FileNotFoundError(f"oracle*.sdf 를 찾을 수 없음: {docking_dir}")
    print(f"[1/3] SDF {len(sdfs)}개에서 후보 추출 중...")
    best: dict[str, float] = {}
    for sdf in sdfs:
        for mol in Chem.SDMolSupplier(sdf):
            if mol is None:
                continue
            try:
                score = float(mol.GetProp("docking_score"))
            except KeyError:
                continue
            smi = Chem.MolToSmiles(mol)
            if smi not in best or score < best[smi]:
                best[smi] = score  # 더 좋은(낮은) 점수 유지
    print(f"      고유 분자 {len(best)}개 수집")
    return best


# ---------- 다양성 기반 top-K (RxnFlow 논문 방식, Tanimoto 0.5) ----------
def diverse_top_k(rows: list[tuple[str, float]], k: int, thresh: float = 0.5):
    rows = sorted(rows, key=lambda r: r[1])  # 점수 낮을수록 좋음
    kept = [rows[0]]
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(rows[0][0]))]
    seen = {rows[0][0]}
    for smi, sc in rows[1:]:
        if smi in seen:
            continue
        seen.add(smi)
        fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smi))
        if max(DataStructs.BulkTanimotoSimilarity(fp, fps)) >= (1 - thresh):
            continue
        fps.append(fp)
        kept.append((smi, sc))
        if len(kept) >= k:
            break
    return kept


# ---------- 2) 재도킹 (verbose 면 배치로 쪼개 진행률 출력) ----------
def redock(vina: VinaReward, smiles_list: list[str], verbose: bool, batch: int) -> list[float]:
    n = len(smiles_list)
    if not verbose:
        return vina.run_smiles(smiles_list)

    scores: list[float] = []
    t0 = time.time()
    for i in range(0, n, batch):
        chunk = smiles_list[i : i + batch]
        scores.extend(vina.run_smiles(chunk))
        done = min(i + batch, n)
        el = time.time() - t0
        eta = (el / done) * (n - done) if done else 0.0
        print(
            f"      재도킹 {done}/{n} ({done / n * 100:.1f}%) | "
            f"경과 {el / 60:.1f}분 | ETA {eta / 60:.1f}분",
            flush=True,
        )
    return scores


def main():
    p = argparse.ArgumentParser(description="이중 타겟 스크리닝 (1차 SDF -> 2차 재도킹 -> 합산)")
    # 1차 타겟 (이미 돌린 run)
    p.add_argument("--docking_dir", required=True, help="1차 타겟 run 의 docking 폴더 (oracle*.sdf 포함)")
    p.add_argument("--target1_name", default="target1", help="1차 타겟 이름 (컬럼 라벨용)")
    # 2차 타겟 (재도킹 대상) -- 하드코딩 없음, 전부 인자
    p.add_argument("--target2_name", default="target2", help="2차 타겟 이름 (컬럼 라벨용)")
    p.add_argument("--target2_protein", required=True, help="2차 타겟 단백질 PDB 경로")
    p.add_argument("--target2_ref_ligand", help="2차 타겟 공결정 리간드 (sdf/mol2/pdb) - 포켓 자동 지정")
    p.add_argument("--target2_center", nargs=3, type=float, help="2차 타겟 포켓 중심좌표 (리간드 대신)")
    p.add_argument("--target2_size", nargs=3, type=float, default=[22.5, 22.5, 22.5], help="도킹 박스 크기")
    p.add_argument("--search_mode", default="fast", choices=["fast", "balance", "detail"])
    # 추림/출력
    p.add_argument("--pre_topk", type=int, default=2000, help="재도킹 전 1차 후보를 몇 개로 추릴지 (0=전부)")
    p.add_argument("--diverse", action="store_true", help="재도킹 전 Tanimoto 0.5 다양성 필터 적용")
    p.add_argument("--final_topk", type=int, default=1000, help="최종 dual hit 개수")
    p.add_argument("--num_workers", type=int, default=None, help="도킹 병렬 워커 수 (기본=가용 코어수)")
    p.add_argument("--out", required=True, help="결과 CSV 경로")
    # verbose
    p.add_argument("--verbose", action="store_true", help="재도킹 진행률/ETA 출력")
    p.add_argument("--batch", type=int, default=500, help="verbose 진행률 출력 배치 크기")
    args = p.parse_args()

    if args.target2_ref_ligand is None and args.target2_center is None:
        p.error("--target2_ref_ligand 또는 --target2_center 중 하나는 필수입니다.")

    t1, t2 = args.target1_name, args.target2_name

    # 1) 추출
    best = extract_candidates(args.docking_dir)
    rows = list(best.items())

    # 재도킹 전 추림 (부담/시간 절약)
    if args.pre_topk > 0:
        rows = diverse_top_k(rows, args.pre_topk) if args.diverse else sorted(rows, key=lambda r: r[1])[: args.pre_topk]
    smiles_list = [smi for smi, _ in rows]
    t1_scores = [sc for _, sc in rows]
    print(f"[2/3] 재도킹 대상 {len(smiles_list)}개 -> {t2} ({args.target2_protein})")

    # 2) 2차 타겟 재도킹 (CLI 로 받은 단백질/포켓)
    vina = VinaReward(
        protein_pdb_path=args.target2_protein,
        center=tuple(args.target2_center) if args.target2_center else None,
        ref_ligand_path=args.target2_ref_ligand,
        size=tuple(args.target2_size),
        search_mode=args.search_mode,
        num_workers=args.num_workers,
    )
    t2_scores = redock(vina, smiles_list, args.verbose, args.batch)

    # 3) 합산 + 정렬
    print("[3/3] 점수 합산 및 정렬...")
    df = pd.DataFrame({
        "smiles": smiles_list,
        f"vina_{t1}": t1_scores,
        f"vina_{t2}": t2_scores,
    })
    df["dual_score"] = df[f"vina_{t1}"] + df[f"vina_{t2}"]  # 낮을수록 좋음
    df = df.sort_values("dual_score").head(args.final_topk).reset_index(drop=True)
    df.insert(0, "hit_id", [f"HIT_{i + 1:04d}" for i in range(len(df))])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n저장 완료: {args.out} ({len(df)} dual-target hits)")
    print(f"  best {t1}: {df[f'vina_{t1}'].min():.2f} | best {t2}: {df[f'vina_{t2}'].min():.2f}")
    print(f"  best dual_score: {df['dual_score'].min():.2f}")


if __name__ == "__main__":
    main()