"""
Module for Drug Recommendation and Report Generation.

This module queries the ChEMBL database to retrieve detailed information about
candidate drugs (clinical phase, structure, etc.) and finds similar molecules
based on SMILES structures. It then generates a formatted text report.
"""
import os
import requests
import json
from urllib.parse import quote

def _chembl_from_pert_name(name):
    """
    Searches ChEMBL for a molecule by its name.

    Args:
        name (str): Drug name to search.

    Returns:
        dict or None: Dictionary containing ChEMBL ID, score, phase, structure, etc.,
                      or None if not found.
    """
    print(f"Searching ChEMBL for: {name}")
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"
    try:
        r = requests.get(url, params={"q": name})
        r.raise_for_status()
        data = r.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data for {name}: {e}")
        return None

    molecules = data.get("molecules", [])
    if not molecules:
        return None
    m0 = molecules[0]
    return {
        "chembl_id": m0["molecule_chembl_id"],
        "score": m0.get("score"),
        # 여기부터는 약물과 관련된 정보
        "atc_classifications": list(set(m0["atc_classifications"])),   # 실제 판매되는 약물의 이름
        "pref_name": m0.get("pref_name"),   # 공식적인 약물 이름
        "max_phase": m0.get("max_phase"),   # 임상시험 단계
        "therapeutic_flag": m0.get("therapeutic_flag"),   # 치료용 약물 여부
        # 여기서부터는 구조
        "smiles": m0['molecule_structures'].get("canonical_smiles") if m0.get('molecule_structures') else None,
        "standard_inchi_key": m0['molecule_structures'].get("standard_inchi_key") if m0.get('molecule_structures') else None,
    }

def _parsing_to_json(data):
    """
    Parses raw ChEMBL molecule data into a simplified dictionary.

    Args:
        data (dict): Raw molecule data from ChEMBL API.

    Returns:
        dict: Simplified dictionary with relevant fields.
    """
    return {
        "chembl_id": data["molecule_chembl_id"],
        "similarity": data.get("similarity"),
        # 여기부터는 약물과 관련된 정보
        "atc_classifications": list(set(data["atc_classifications"])),   # 실제 판매되는 약물의 이름, 없으면 판매중 아님
        "pref_name": data.get("pref_name"),   # 공식적인 약물 이름
        "max_phase": data.get("max_phase"),   # 임상시험 단계
        "therapeutic_flag": data.get("therapeutic_flag"),   # 치료용 약물 여부
        # 이거는 구조 기반 검색을 위함
        "smiles": data['molecule_structures'].get("canonical_smiles") if data.get('molecule_structures') else None,
        "standard_inchi_key": data['molecule_structures'].get("standard_inchi_key") if data.get('molecule_structures') else None,
        }

def _similar_from_smile(smile):
    """
    Searches for molecules with similar structures using SMILES.

    Args:
        smile (str): Canonical SMILES string.

    Returns:
        list or None: List of similar molecules (parsed), or None if error/not found.
    """
    print(f"Searching similar molecules for SMILES: {smile}")
    score = 60
    ls_num = 20

    # SMILES에 특수문자(# 등)가 포함될 수 있으므로 URL 인코딩 필요
    encoded_smile = quote(smile)
    url  = f"https://www.ebi.ac.uk/chembl/api/data/similarity/{encoded_smile}/{score}?format=json&limit={ls_num}"

    try:
        r = requests.get(url)
        r.raise_for_status()
        data = r.json()
        return [_parsing_to_json(item) for item in data['molecules']]
    except requests.exceptions.RequestException as e:
        print(f"Error fetching similar molecules for SMILES {smile}: {e}")
        return None

def _save_to_file(results, new_brd, out_path):
    """
    Saves the recommendation results to a formatted text file.

    Writes an abstract section with a tree-like structure of similar drugs
    and a detailed information section for each candidate.

    Args:
        results (dict): Dictionary of ChEMBL search results.
        new_brd (dict): Dictionary of original candidate drugs.
        out_path (str): Output file path.
    """
    # 결과 파일로 저장
    symbol = ['├──','└──']
    with open(out_path, 'w', encoding='utf-8') as f:
        # 개괄
        f.write("<Abstract>\n")
        for idx, (key, val) in enumerate(results.items()):
            info = new_brd[key]
            f.write(f"{idx+1}. {info['name']}({val['chembl_id']}) [TAG: {info['tag']:.2f}]\n")
            ls = []
            if val['similar_molecules'] is None:
                f.write("\n")
                continue
            for similar in val["similar_molecules"]:
                if similar['chembl_id'] == val['chembl_id']:
                    continue
                if similar['max_phase'] is None:
                    continue
                if float(similar['max_phase']) < 2:
                    continue
                ls.append(f"{similar['pref_name']}({similar['max_phase']})")
            for idx, chembl in enumerate(ls):
                if idx == len(ls) - 1:
                    f.write(f"{symbol[1]}{chembl}\n")
                else:
                    f.write(f"{symbol[0]}{chembl}\n")
            f.write("\n")
                
        
        # 세부 정보 작성
        f.write("\n\n<Detailed Information>\n")
        for idx, (key, val) in enumerate(results.items()):
            info = new_brd[key]
            f.write(f"===후보 물질 {idx+1}: {info['name']}===\n")
            f.write(f"TAG Score: {info['tag']:.2f}\n")
            f.write(f"ChEMBL Accession ID: {val['chembl_id']}\n")
            f.write(f"SMILES: {val['smiles']}\n")
            f.write(f"Standard InChI Key: {val['standard_inchi_key']}\n")
            if val['max_phase'] is None:
                f.write("임상 시험 정보가 없습니다.")
                f.write("\n\n\n")
                continue
            f.write("약물 정보\n")
            f.write(f"- 공식 약물 이름: {val['pref_name']}\n")
            if val['atc_classifications']:
                f.write(f"- 실제 판매되는 이름: {', '.join(val['atc_classifications'])}\n")
            f.write(f"- 최대 임상시험 단계: {val['max_phase']}\n")
            if val['therapeutic_flag']:
                f.write("- 실제 치료에 사용되고 있습니다.\n")
            else:
                if val['max_phase'] == 0:
                    f.write("- 실제 치료에 사용되지 않습니다.\n")
                else:
                    f.write("- 아직 임상단계로 상용화되지 않았습니다.\n")

            if val['similar_molecules'] is None:
                f.write("\n\n")
                continue
            idx = 1
            for similar in val["similar_molecules"]:
                if similar['chembl_id'] == val['chembl_id']:
                    continue
                if similar['max_phase'] is None:
                    continue
                if float(similar['max_phase']) < 2:
                    continue
                f.write(f"\n\t후보 물질 {key}과 유사물질 {idx}번\n")
                f.write(f"\t후보 물질과의 유사도 점수: {float(similar['similarity']):.1f}%\n")
                f.write(f"\tChEMBL Accession ID: {similar['chembl_id']}\n")
                f.write(f"\tSMILES: {similar['smiles']}\n")
                f.write(f"\tStandard InChI Key: {similar['standard_inchi_key']}\n")
                f.write("\t약물 정보\n")
                f.write(f"\t- 공식 약물 이름: {similar['pref_name']}\n")
                if similar['atc_classifications']:
                    f.write(f"\t- 실제 판매되는 이름: {', '.join(similar['atc_classifications'])}\n")
                f.write(f"\t- 최대 임상시험 단계: {similar['max_phase']}\n")
                if similar['therapeutic_flag']:
                    f.write("\t- 실제 치료에 사용되고 있습니다.\n")
                else:
                    if similar['max_phase'] == "None":
                        f.write("\t- 치료에 사용될 수 없습니다.\n")
                    else:
                        f.write("\t- 아직 임상단계로 상용화되지 않았습니다.\n")
                idx += 1
            f.write("\n\n")

def make_result(brd_dict, work_path):
    """
    Main function to generate the recommendation report.

    Iterates through candidate drugs, queries ChEMBL for details and similar molecules,
    and saves the consolidated results to a file.

    Args:
        brd_dict (dict): Dictionary of candidate drugs from CMap analysis.
        work_path (str): Output file path.
    """
    results = {}
    new_brd = brd_dict.copy()

    # brd id로 ChEMBL 검색
    for brd, info in brd_dict.items():
        name = info['name']
        temp = _chembl_from_pert_name(name)
        if temp is None:
            continue
        # 중복 제거(다른 이름인데 동일한 ChEMBL_id인 경우)
        if temp["chembl_id"] in [res['chembl_id'] for res in results.values()]:
            new_brd.pop(brd)
            continue
        results[brd] = temp

    # Similar molecules 검색
    for key, _ in new_brd.items():
        smile = results[key]['smiles']
        if smile:
            results[key]['similar_molecules'] = _similar_from_smile(smile)
        else:
            results[key]['similar_molecules'] = None
    
    # 결과 파일로 저장
    _save_to_file(results, new_brd, work_path)
