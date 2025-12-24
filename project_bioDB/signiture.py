"""
Module for Gene Signature Extraction.

This module handles the reading of gene expression data, grouping of samples (WT vs MUT),
normalization using PyDESeq2, and extraction of differentially expressed genes (DEGs)
to form the gene signature (Up/Down regulated genes).
"""
import pandas as pd
import os
from pydeseq2.default_inference import DefaultInference
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def _read_col_file(home):
    """
    Reads the dataset label file to identify sample groups.

    Args:
        home (str): Home directory path containing 'dataset_label.txt'.

    Returns:
        dict: A dictionary mapping filenames to lists of group labels (e.g., 'mut', 'wt').
    """
    file_groups = {}
    target = os.path.join(home, 'dataset_label.txt')
    with open(target, 'r', encoding='utf-8') as f:
        current_key = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                continue
            elif line.startswith('>>'):
                current_key = line.replace('>>', '').strip()
            else:
                if current_key:
                    value = line.split(',')
                    # 대소문자 통일 (mut, wt) 및 공백 제거
                    file_groups[current_key] = [val.strip().lower() for val in value]
    return file_groups  
            

def _regularization_protocol(path, file_groups):
    """
    Performs normalization and DEG analysis for multiple files combined.

    Used when the '--same' protocol is selected. Merges dataframes, filters low counts,
    and runs DESeq2 analysis.

    Args:
        path (str): Directory path containing expression files.
        file_groups (dict): Dictionary of filenames and their group labels.

    Returns:
        pd.DataFrame: DESeq2 results dataframe.
    """
    dfs = []
    group = []
    experiment = []

    for fname, groups in file_groups.items():
        file_path = os.path.join(path, fname)
        if fname.endswith('.tsv'):
            df = pd.read_csv(file_path, sep="\t", index_col=0)
        elif fname.endswith('.csv'):
            df = pd.read_csv(file_path, sep=",", index_col=0)
        
        dfs.append(df)
        group.extend(groups)
        experiment.extend([fname] * df.shape[1])

    counts_df = pd.concat(dfs, axis=1)
    counts_df = counts_df.fillna(0) #na는 0으로
    counts_df = counts_df.T
    counts_df = counts_df[~(counts_df.lt(10).all(axis=1))] #10개 이하는 지워준다.

    meta_data = pd.DataFrame({
        "group": group,
        "experiment": experiment
    }, index=counts_df.index)

    inference = DefaultInference(n_cpus=1)

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=meta_data,
        design="~ experiment + group",
        refit_cooks=True,
        inference=inference
    )

    dds.deseq2()

    ds = DeseqStats(
        dds,
        contrast=["group", "mut", "wt"],
        inference=inference
    )

    ds.summary()
    res = ds.results_df
    return res



def _signiture_analysis(res):
    """
    Filters DESeq2 results to extract top Up and Down regulated genes.

    Args:
        res (pd.DataFrame): DESeq2 results dataframe.

    Returns:
        tuple: Two lists containing Up-regulated and Down-regulated gene symbols.
    """
    N = 30   # CMap에 넣을 up / down 유전자 개수
    
    # 1) 유의한 DEG만 보고 싶으면 padj 기준 0.1 이하 필터링
    deg = res[res["padj"] < 0.1]
    print(f"유의한 DEG 수: {deg.shape[0]}")

    # 2) CMap용 Up / Down 상위 N개 선택
    up = deg[deg["log2FoldChange"] > 0].sort_values("padj", ascending=True)
    down = deg[deg["log2FoldChange"] < 0].sort_values("padj", ascending=True)

    # 원하는 갯수 미달일 경우 padj에서 상위 N개로 채우기
    if len(up) < N:
        up = res[res["log2FoldChange"] > 0].sort_values("padj", ascending=True)[:N]
    else:
        up = up[:N]  # 너무 많으면 N개로 자름
        
    if len(down) < N:
        down = res[res["log2FoldChange"] < 0].sort_values("padj", ascending=True)[:N]
    else:
        down = down[:N] # 너무 많으면 N개로 자름

    # pandas -> list
    up_genes = up.index.tolist()
    down_genes = down.index.tolist()

    print(f"Final Up genes: {len(up_genes)}, Down genes: {len(down_genes)}")
    
    return up_genes, down_genes



def _regularization_single_protocol(path, file_groups):
    """
    Performs normalization and DEG analysis for each file individually.

    Used when the '--same' protocol is NOT selected. Runs DESeq2 for each file
    and collects significant genes.

    Args:
        path (str): Directory path containing expression files.
        file_groups (dict): Dictionary of filenames and their group labels.

    Returns:
        tuple: Lists of all Up and Down genes found across files.
    """
    up_gene = []
    down_gene = []
    for fname, groups in file_groups.items():
        file_path = os.path.join(path, fname)
        if fname.endswith('.tsv'):
            df = pd.read_csv(file_path, sep="\t", index_col=0)
        elif fname.endswith('.csv'):
            df = pd.read_csv(file_path, sep=",", index_col=0)

        counts_df = df.fillna(0) #na는 0으로
        counts_df = counts_df.T
        counts_df = counts_df[~(counts_df.lt(10).all(axis=1))] #10개 이하는 지워준다.

        meta_data = pd.DataFrame({
            "group": groups,
            "experiment": [fname] * df.shape[1]
        }, index=counts_df.index)

        inference = DefaultInference(n_cpus=1)
        
        dds = DeseqDataSet(
            counts=counts_df,
            metadata=meta_data,
            design="~ experiment + group",
            refit_cooks=True,
            inference=inference
        )

        dds.deseq2()
        ds = DeseqStats(
            dds,
            contrast=["group", "mut", "wt"],
            inference=inference
        )

        ds.summary()
        res = ds.results_df
        sig_up, sig_down = _signiture_analysis(res)
        up_gene.extend(sig_up)
        down_gene.extend(sig_down)
    return up_gene, down_gene



def extract_signiture(same, path):
    """
    Main function to extract gene signatures.

    Orchestrates the reading of labels and execution of either the combined
    or single regularization protocol based on the 'same' flag.

    Args:
        same (bool): Flag indicating if the same protocol applies to all files.
        path (str): Working directory path.

    Returns:
        tuple: Final lists of Up and Down regulated gene symbols.
    """
    file_groups = _read_col_file(path)

    if same:  # 동일한 실험 프로토콜 -> 한번에 합쳐서 진행
        res = _regularization_protocol(path, file_groups)
        up_genes, down_genes = _signiture_analysis(res)
    else:     # 전혀 다른 프로토콜 -> 각각 분석 후 합침
        up_gene, down_gene = _regularization_single_protocol(path, file_groups)

        up_genes = [g for g in up_gene if up_gene.count(g) >= len(file_groups) / 2]
        down_genes = [g for g in down_gene if down_gene.count(g) >= len(file_groups) / 2]
    
    return up_genes, down_genes
