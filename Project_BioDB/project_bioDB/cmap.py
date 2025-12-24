"""
Module for CMap API Interaction.

This module handles interactions with the CMap (Connectivity Map) API, including
mapping gene symbols to Entrez IDs, submitting analysis jobs, monitoring status,
and downloading results.
"""
import requests
import gzip
import os
import pandas as pd
import time

def _get_L1000_BING_genes(path):
    """
    Downloads and parses the L1000 BING gene list.

    Args:
        path (str): Directory to save/load the gene info file.

    Returns:
        dict: Dictionary mapping gene symbols to Entrez IDs.
    """
    out_fname = 'L1000_BING_genes.txt.gz'
    out_path = os.path.join(path, out_fname)

    if not os.path.exists(out_path):
        print(f"Downloading gene_info from GEO...")
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_gene_info.txt.gz"

        r = requests.get(url)
        r.raise_for_status()

        # 1) gzip 파일로 저장
        with open(out_path, "wb") as f:
            f.write(r.content)
        print("Download complete:", out_path)

    # 2) gzip 압축 풀기
    with gzip.open(out_path, 'rt') as f_in:
        gene_info = pd.read_csv(f_in, sep='\t')

    # 3) BING 유전자만 필터링
    bing = gene_info[gene_info['pr_is_bing'] == 1]
    data = bing.iloc[:, 0:2]  # 'pr_gene_id', 'pr_gene_symbol' 열 선택

    entrez_ids = {}
    for entrez, gene_symbol in data.values:
        entrez_ids[gene_symbol] = str(entrez)
    return entrez_ids

def _find_another_symbol(gene_symbol):
    """
    Searches for alternative gene symbols using the genenames.org API.

    Args:
        gene_symbol (str): Gene symbol to search.

    Returns:
        list: List of previous or alias symbols found.
    """
    # fetch/symbol을 사용하면 해당 유전자의 모든 정보(prev, alias 포함)를 가져옵니다.
    headers = {"Accept": "application/json"}
    
    r = requests.get(f"https://rest.genenames.org/fetch/symbol/{gene_symbol}", headers=headers)
    
    alt_symbol = []
    if r.ok:
        data = r.json()
        if data['response']['numFound'] > 0:
            doc = data['response']['docs'][0]
            alt_symbol.extend(doc.get('prev_symbol', []))
            alt_symbol.extend(doc.get('alias_symbol', []))
    return alt_symbol

def _transfer_to_entrez(up_genes, down_genes, entrez_ids):
    """
    Converts gene symbols to Entrez IDs.

    Tries to map directly, and falls back to searching for alternative symbols
    if a direct match is not found. Logs the matching process to 'gene_match_log.txt'.

    Args:
        up_genes (list): List of Up-regulated gene symbols.
        down_genes (list): List of Down-regulated gene symbols.
        entrez_ids (dict): Dictionary of valid L1000 Entrez IDs.

    Returns:
        tuple: Lists of mapped Entrez IDs for Up and Down genes.
    """
    up_entrez = []
    for up in up_genes:
        if up in entrez_ids:
            up_entrez.append(entrez_ids[up])
        else:
            alt_symbols = _find_another_symbol(up)
            for alt_symbol in alt_symbols:
                if alt_symbol and alt_symbol in entrez_ids:
                    up_entrez.append(entrez_ids[alt_symbol])
                    break

    down_entrez = []
    for down in down_genes:
        if down in entrez_ids:
            down_entrez.append(entrez_ids[down])
        else:
            alt_symbols = _find_another_symbol(down)
            for alt_symbol in alt_symbols:
                if alt_symbol and alt_symbol in entrez_ids:
                    down_entrez.append(entrez_ids[alt_symbol])
                    break
    
    print(f"Entrez ID 매칭 성공: Up: {len(up_entrez)}, Down: {len(down_entrez)}")
    return up_entrez, down_entrez

def _request_cmap(up_entrez, down_entrez, api_key):
    """
    Submits a job to the CMap API.

    Args:
        up_entrez (list): List of Up-regulated Entrez IDs.
        down_entrez (list): List of Down-regulated Entrez IDs.
        api_key (str): CMap user API key.

    Returns:
        str or None: The Job ID if successful, None otherwise.
    """
    print(f"Entrez ID 매칭 성공: Up: {len(up_entrez)}, Down: {len(down_entrez)}")
    up_line = "TAG\t\t" + "\t".join(up_entrez)
    dn_line = "TAG\t\t" + "\t".join(down_entrez)

    url = 'https://api.clue.io/api/jobs'

    headers = {
        'user_key': api_key,
        'Content-Type': 'application/json',
        'Accept': 'application/json'
    }
    payload = {
        "tool_id": "sig_gutc_tool",
        "name": "Sample_MUT_vs_WT_test",
        "data_type": "L1000",
        "dataset": "Touchstone",
        "ignoreWarnings": True,          # 경고 무시하고 진행
        "uptag-cmapfile": up_line,
        "dntag-cmapfile": dn_line,
    }
    response = requests.post(url, headers=headers, json=payload)
    print(f"Status Code: {response.status_code}")

    job_id = None

    if response.status_code in (200, 201, 202):
        data = response.json()
        print("Job meta:", data)
        job_id = data.get("result", {}).get("job_id")
        print("Job ID:", job_id)

    else:
        print("ERROR")
        print("status:", response.status_code)
    return job_id

def _status(job_id, api_key):
    """
    Checks the status of a CMap job.

    Args:
        job_id (str): The Job ID to check.
        api_key (str): CMap user API key.

    Returns:
        bool: True if completed, False if submitted/pending.
    """
    status_url = f"https://api.clue.io/api/jobs/findByJobId/{job_id}"
    status_resp = requests.get(status_url, headers={"user_key": api_key})
    r = status_resp.json()
    status = r.get('status')
    if status == 'submitted' or status == 'pending':
        print(f"{status}... 3분 대기")
        return False
    elif status == 'completed':
        return True


def _get_cmap_result(out_path, job_id, api_key):
    """
    Monitors job status and downloads the result when complete.

    Waits in 3-minute intervals until the job is done, then downloads
    the result tar.gz file.

    Args:
        out_path (str): Directory to save the result.
        job_id (str): CMap Job ID.
        api_key (str): CMap user API key.

    Returns:
        str: The filename of the downloaded result.
    """
    print(f"Job {job_id} 상태 모니터링 시작 (3분 간격)...")
    while True:
        if _status(job_id, api_key):
            break
        time.sleep(180)  # 3분 대기

    url = f"https://api.clue.io/api/jobs/findByJobId/{job_id}"
    resp = requests.get(url, headers={"user_key": api_key})
    print("Status code:", resp.status_code)
    print("Status body:", resp.text)

    standard_result = resp.json().get("download_url")
    print("raw download_url field:", standard_result)

    if standard_result:
        # 앞에 프로토콜 붙이기
        if standard_result.startswith("//"):
            download_url = "https:" + standard_result
        else:
            download_url = standard_result

        print("Download URL:", download_url)

        r = requests.get(download_url, stream=True)
        r.raise_for_status() 
        
        out_fname = f"cmap_result_{job_id}.tar.gz"
        result_path = os.path.join(out_path, out_fname)

        with open(result_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        print("CMap result downloaded:", out_path)
    else:
        print("Download_url 필드가 없음")
    return out_fname

def cmap_analysis(api_key, up_genes, down_genes, work_path):
    """
    Main function to orchestrate CMap analysis.

    Gets L1000 genes, maps input genes to Entrez IDs, submits the job,
    and downloads the result.

    Args:
        api_key (str): CMap API key.
        up_genes (list): Up-regulated gene symbols.
        down_genes (list): Down-regulated gene symbols.
        work_path (str): Working directory.

    Returns:
        str: The filename of the downloaded CMap result.
    """
    entrez_ids = _get_L1000_BING_genes(work_path)
    up_entrez, down_entrez = _transfer_to_entrez(up_genes, down_genes, entrez_ids)
    job_id = _request_cmap(up_entrez, down_entrez, api_key)
    
    if not job_id:
        raise ValueError("CMap API 요청 실패: Job ID를 받아오지 못했습니다. 위 에러 로그를 확인하세요.")

    cmap_name = _get_cmap_result(work_path, job_id, api_key)
    return cmap_name
