"""
Main execution module for Project_BioDB.

This module serves as the entry point for the pipeline. It handles argument parsing,
directory setup, and orchestrates the workflow from gene signature extraction to
drug recommendation generation.
"""
from project_bioDB import *
import argparse
import pandas as pd
import os

from project_bioDB.cmap import _get_L1000_BING_genes, _transfer_to_entrez

def main():
    """
    Main function to execute the Project_BioDB pipeline.

    Parses command-line arguments, sets up the working directory, validates input files,
    and sequentially calls the analysis modules:
    1. Extract gene signatures (Up/Down regulated genes).
    2. Map genes to Entrez IDs.
    3. Run CMap analysis (API request and result download).
    4. Process CMap results to find candidate drugs.
    5. Generate final recommendations with ChEMBL data.
    """
    # 현재 파일(main.py)이 있는 디렉토리를 작업 디렉토리로 설정
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print(f"Current working directory: {os.getcwd()}")

    parser = argparse.ArgumentParser(description="Process gene expression files.")
    # dataset_label.txt를 사용하므로 files 인자 제거
    parser.add_argument('--same', action='store_true', help="Use same protocol for all files (default: False)")
    parser.add_argument('--api_key', type=str, required=True, help="CMap API Key")
    
    args = parser.parse_args()

    home = os.getcwd()
    work_path = os.path.join(home, 'data')
    
    # data 폴더 생성
    os.makedirs(work_path, exist_ok=True)

    # dataset_label.txt 유효성 검사
    target = os.path.join(home, 'dataset_label.txt')
    with open(target, 'r', encoding='utf-8') as f: 
        for line in f:
            if line.startswith('#'):
                continue
            elif line.startswith('>>'):
                break
            raise ValueError("dataset_label.txt 파일에 '>>'로 시작하는 라인이 없습니다.")
        

    # python main.py (--same) --api_key YOUR_KEY
    print(f"Protocol: {'Same' if args.same else 'Different'}")
    
    # 1) Signature 추출
    up_genes, down_genes = extract_signiture(args.same, home)
    print(f"Extracted Signature: {len(up_genes)} Up, {len(down_genes)} Down")
    
    # 2) CMap 분석
    result_name = cmap_analysis(args.api_key, up_genes, down_genes, work_path)
    
    # 3) CMap 결과 처리
    final_compound = get_drug_list(work_path, result_name)

    # 4) 결과 추천
    output_fname = result_name.replace('cmap_result', 'Recommendations').replace('.tar.gz', '.txt')
    output_path = os.path.join(work_path, output_fname)
    make_result(final_compound, output_path)
    print(f"Analysis saved to {output_path}")
    

if __name__ == "__main__":
    main()
