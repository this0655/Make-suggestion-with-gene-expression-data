"""
Module for analyzing CMap results.

This module handles the decompression of CMap result files, parsing of GCT files
(Gene Cluster Text), and extraction of top drug candidates based on Connectivity Scores (TAG).
"""
import tarfile
import os
from cmapPy.pandasGEXpress.parse_gct import parse

def _decompression_cmap_results(work_path, result_name):
    """
    Decompresses the CMap result tar.gz file.

    Args:
        work_path (str): The directory where the result file is located.
        result_name (str): The name of the tar.gz file.

    Returns:
        str: The name of the top-level directory extracted from the archive.
    """
    # 1. tar.gz 압축 해제
    result_path = os.path.join(work_path, result_name)
    with tarfile.open(result_path, 'r:gz') as tar:
        all_members = tar.getmembers()
        first_path = all_members[0].name if all_members else 'No members found'
        tar.extractall(work_path)
    return first_path



def _read_gct(base_dir, fname):
    """
    Reads a GCT file and returns its data and metadata.

    Args:
        base_dir (str): Base directory path.
        fname (str): GCT filename.

    Returns:
        dict or None: A dictionary containing 'df' (data), 'row_metadata', and 'col_metadata',
                      or None if the file does not exist.
    """
    result = {}
    path = os.path.join(base_dir, fname)
    if os.path.exists(path):
        print(f"Reading from {path}")
        gct = parse(path)
        result['df'] = gct.data_df  # DataFrame
        result["row_metadata"] = gct.row_metadata_df
        result["col_metadata"] = gct.col_metadata_df
        return result
    else:
        print(f"File not found: {path}")
        return None



def _get_pert_summary(work_path, first_path):
    """
    Retrieves the perturbation summary GCT file.

    Args:
        work_path (str): Working directory path.
        first_path (str): Extracted directory name.

    Returns:
        dict: Parsed GCT data for pert_id_summary.gct.
    """
    # arfs/TAG/pert_id_summary.gct 경로 구성
    target_dir = os.path.join(work_path, first_path, 'arfs', 'TAG')
    return _read_gct(target_dir, 'pert_id_summary.gct')



def _get_csn(work_path, first_path):
    """
    Searches for and retrieves the connectivity score (cs) GCT file.

    Args:
        work_path (str): Working directory path.
        first_path (str): Extracted directory name.

    Returns:
        dict or None: Parsed GCT data for the found cs_*.gct file.
    """
    # cs_*.gct 파일 찾기 (동적 탐색)
    base_dir = os.path.join(work_path, first_path)
    ls = os.listdir(base_dir)
    for file in ls:
        if file.startswith('cs_') and file.endswith('.gct'):
            return _read_gct(base_dir, file)
    return None



def get_drug_list(work_path, result_name):
    """
    Extracts the top drug candidates from CMap results.

    Decompresses the result, reads the summary and connectivity score files,
    filters for Broad ID (BRD) compounds, and selects the top 10 drugs with
    the lowest TAG scores (indicating reversal of disease signature).

    Args:
        work_path (str): Working directory path.
        result_name (str): CMap result filename.

    Returns:
        dict: A dictionary mapping Broad IDs to drug names and TAG scores.
    """
    # 1. 결과 압축 해제 및 첫 디렉토리 경로 얻기
    first_path = _decompression_cmap_results(work_path, result_name)
    
    # 2. pert_id_summary.gct 파일 읽기
    pert_summary = _get_pert_summary(work_path, first_path)
    
    # 3. cs_*.gct 파일 읽기 (동적 탐색)
    csn = _get_csn(work_path, first_path)

    df1 = pert_summary['df']
    df1 = df1[df1.index.str.startswith('BRD')]
    # TAG 점수가 낮은 순서대로 상위 10개 추출
    df1 = df1.sort_values('TAG', ascending=True).head(10)
    drug_list = df1.index.to_list()

    df2 = csn['row_metadata']
    drug_dict = {}
    for brd in drug_list:
        row = df2['pert_iname'].loc[df2['pert_id'] == brd]
        tag_score = df1.loc[brd, 'TAG']
        drug_dict[brd] = {'name': row.values[0], 'tag': float(tag_score)}
    return drug_dict

