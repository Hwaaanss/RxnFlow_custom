# 0) Git Clone
```bash
git clone https://github.com/SeonghwanSeo/RxnFlow.git
cd RxnFlow
```


# 1) Setting Environment
```bash
conda new -s rxnflow python=3.12
conda activate rxnflow

pip install -e '.[unidock,pmnet,dev]' --find-links https://data.pyg.org/whl/torch-2.5.1+cu121.html
conda install -c conda-forge unidock==1.1.2 

pip install gdown
gdown --folder https://drive.google.com/drive/folders/1e5pPZaTRGhvEMky3K2OKQ9-jV_NweK-a
gdown --folder https://drive.google.com/drive/folders/1UDAxFPYWWOl_GJK_DMduoFLPP-zMEf0f
```

zincfrag*.smi.gz 파일들은 ./data/building_blocks/ 로 이동
LIT-PCBA.tar.gz 과 CrossDocked2020_*.tar.gz 파일들은 ./data/experiments/ 로 이동


# 2) Data Refine
```bash
cd data
tar -xzf ./experiments/LIT-PCBA.tar.gz
tar -xzf ./experiments/CrossDocked2020_all.tar.gz
tar -xzf ./experiments/CrossDocked2020_test.tar.gz

zcat ./building_blocks/zincfrag_5.6M.smi.gz > ./building_blocks/zincfrag_5.6M.smi

python scripts/a_refine_smi.py \
  -b building_blocks/zincfrag_5.6M.smi \
  -o building_blocks/zincfrag_5.6M_clean.smi \
  --cpu 8

python scripts/b_create_env.py \
  -b building_blocks/zincfrag_5.6M_clean.smi \
  -o ./envs/zincfrag-5p6m-clean \
  -t ./templates/real.txt --cpu 8
cd ..
```

-o 뒤 경로는 신규 생성 디렉토리이고, exist_ok=False 라서 폴더명 중복 시 에러
./rcsb 폴더명을 targets로 변경


# 3) Extract Target Ligand
```bash
echo "=== 6VQN HETATM ===" && grep "^HETATM" targets/6vqn.pdb | awk '{print $4}' | sort -u
echo "=== 1VJY HETATM ===" && grep "^HETATM" targets/1vjy.pdb | awk '{print $4}' | sort -u
echo "=== 6TVG HETATM ===" && grep "^HETATM" targets/6tvg.pdb | awk '{print $4}' | sort -u
echo "=== 6DHB HETATM ===" && grep "^HETATM" targets/6dhb.pdb | awk '{print $4}' | sort -u
```

위 코드의 결과 중 물(HOH)을 제외한 3글자 출력 약물 리간드를 확보 후 아래 코드의 <LIG_*> 부분에 각각 교체

```bash
grep -E "^HETATM.* <LIG_6VQN> " targets/6vqn.pdb > targets/ligand_6vqn.pdb
grep -E "^HETATM.* <LIG_1VJY> " targets/1vjy.pdb > targets/ligand_1vjy.pdb
grep -E "^HETATM.* <LIG_6TVG> " targets/6tvg.pdb > targets/ligand_6tvg.pdb
grep -E "^HETATM.* <LIG_6DHB> " targets/6dhb.pdb > targets/ligand_6dhb.pdb
```


# 4) Generative Screening
## PD-L1 screening
```bash
python scripts/opt_unidock_moo.py \
  --env_dir ./data/envs/zincfrag-5p6m-clean \
  -o ./logs/pdl1_6vqn \
  -p ./targets/6vqn.pdb \
  -l ./targets/ligand_6vqn.pdb \
  --pretrained_model qed-unif-0-64 \
  --num_iterations 1000
```

## check output's column name and directory
```bash
ls ./logs/pdl1_6vqn/
```

## TGF-BR1 re-docking, dual-target integration and sort
```bash
python ./scripts/dual_target_screen.py \
  --docking_dir ./logs/pdl1_6vqn/docking \
  --target1_name PD-L1 --target2_name TGF-bR1 \
  --target2_protein ./targets/1vjy.pdb \
  --target2_ref_ligand ./targets/ligand_1vjy.pdb \
  --out ./logs/dual_target_hits_6vqn_1vjy.csv \
  --pre_topk 0 --diverse --final_topk 1000 \
  --num_workers 6 --verbose --batch 1000
  ```
  
## 6TVG re-docking, dual-target integration and sort
```bash
python ./scripts/dual_target_screen.py \
  --docking_dir ./logs/pdl1_6vqn/docking \
  --target1_name PD-L1 --target2_name 6TVG \
  --target2_protein ./targets/6tvg.pdb \
  --target2_ref_ligand ./targets/ligand_6tvg.pdb \
  --out ./logs/dual_target_hits_6vqn_6tvg.csv \
  --pre_topk 0 --diverse --final_topk 1000 \
  --num_workers 6 --verbose --batch 1000
  ```

  ## 6DHB re-docking, dual-target integration and sort
```bash
python ./scripts/dual_target_screen.py \
  --docking_dir ./logs/pdl1_6vqn/docking \
  --target1_name PD-L1 --target2_name 6DHB \
  --target2_protein ./targets/6dhb.pdb \
  --target2_ref_ligand ./targets/ligand_6dhb.pdb \
  --out ./logs/dual_target_hits_6vqn_6dhb.csv \
  --pre_topk 0 --diverse --final_topk 1000 \
  --num_workers 6 --verbose --batch 1000
  ```