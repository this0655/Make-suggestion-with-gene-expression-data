# Project_BioDB Code Structure

## Directory Structure (Summary)

```text
Project_BioDB
├── Main.py
│   └── main()
└── project_bioDB/
    ├── signiture.py
    │   ├── extract_signiture(same, path)
    │   ├── _read_col_file(home)
    │   ├── _regularization_protocol(path, file_groups)
    │   ├── _regularization_single_protocol(path, file_groups)
    │   └── _signiture_analysis(res)
    ├── cmap.py
    │   ├── cmap_analysis(api_key, up_genes, down_genes, work_path)
    │   ├── _get_L1000_BING_genes(path)
    │   ├── _transfer_to_entrez(up_genes, down_genes, entrez_ids)
    │   ├── _request_cmap(up_entrez, down_entrez, api_key)
    │   ├── _status(job_id, api_key)
    │   └── _get_cmap_result(out_path, job_id, api_key)
    ├── analyze.py
    │   ├── get_drug_list(work_path, result_name)
    │   ├── _decompression_cmap_results(work_path, result_name)
    │   ├── _get_pert_summary(work_path, first_path)
    │   └── _get_csn(work_path, first_path)
    └── recommendation.py
        ├── make_result(brd_dict, work_path)
        ├── _chembl_from_pert_name(name)
        ├── _similar_from_smile(smile)
        └── _save_to_file(results, new_brd, out_path)
```

## Detailed Description

```text
Project_BioDB
├── Main.py
│   └── main()
│       : 프로그램의 진입점(Entry Point).
│       : 인자 파싱, 데이터 폴더 생성, 전체 파이프라인 실행을 담당합니다.
│
└── project_bioDB/
    ├── signiture.py
    │   ├── extract_signiture(same, path)
    │   │   : Gene Expression 데이터에서 Signature(Up/Down Genes)를 추출하는 메인 함수.
    │   ├── _read_col_file(home)
    │   │   : dataset_label.txt를 읽어 파일별 그룹 정보를 파싱합니다.
    │   ├── _regularization_protocol(path, file_groups)
    │   │   : 모든 파일을 통합하여 정규화 및 DEG 분석을 수행합니다.
    │   ├── _regularization_single_protocol(path, file_groups)
    │   │   : 각 파일을 개별적으로 정규화 및 분석합니다.
    │   └── _signiture_analysis(res)
    │       : DEG 분석 결과에서 유의한 유전자(Up/Down)를 필터링하여 추출합니다.
    │
    ├── cmap.py
    │   ├── cmap_analysis(api_key, up_genes, down_genes, work_path)
    │   │   : CMap API 연동을 총괄하는 함수.
    │   ├── _get_L1000_BING_genes(path)
    │   │   : L1000 유전자 정보를 다운로드하고 BING 유전자만 필터링합니다.
    │   ├── _transfer_to_entrez(up_genes, down_genes, entrez_ids)
    │   │   : 유전자 심볼을 Entrez ID로 변환합니다.
    │   ├── _request_cmap(up_entrez, down_entrez, api_key)
    │   │   : CMap API에 분석 작업을 요청(POST)하고 Job ID를 받습니다.
    │   ├── _status(job_id, api_key)
    │   │   : 작업 상태를 확인합니다.
    │   └── _get_cmap_result(out_path, job_id, api_key)
    │       : 작업 완료 시 결과를 다운로드합니다.
    │
    ├── analyze.py
    │   ├── get_drug_list(work_path, result_name)
    │   │   : CMap 결과 파일에서 상위 약물 후보를 추출합니다.
    │   ├── _decompression_cmap_results(work_path, result_name)
    │   │   : .tar.gz 결과 파일의 압축을 해제합니다.
    │   ├── _get_pert_summary(work_path, first_path)
    │   │   : pert_id_summary.gct 파일을 읽습니다.
    │   └── _get_csn(work_path, first_path)
    │       : cs_*.gct 파일을 동적으로 찾아 읽습니다.
    │
    └── recommendation.py
        ├── make_result(brd_dict, work_path)
        │   : 최종 결과 리포트를 생성하는 메인 함수.
        ├── _chembl_from_pert_name(name)
        │   : 약물 이름으로 ChEMBL 정보를 검색합니다.
        ├── _similar_from_smile(smile)
        │   : SMILES 구조를 기반으로 유사 화합물을 검색합니다.
        └── _save_to_file(results, new_brd, out_path)
            : 수집된 정보를 포맷팅하여 텍스트 파일로 저장합니다.
```
