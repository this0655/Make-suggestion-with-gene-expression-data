# CMap Result File Structure

CMap(clue.io) 분석 결과 파일의 상세 구조입니다.

```text
my_analysis.sig_gutc_tool.{job_id}
├──cs_n1x476251.gct
├──query_config.yaml
├──uptag.gmt
├──arfs
│   ├──index.txt
│   └──TAG
│       ├──column_info.txt
│       ├──pcl_cell.gct
│       ├──pcl_summary.gct
│       ├──pert_id_cell.gct
│       ├──pert_id_summary.gct
│       └──query_info.txt
└──matrices
    ├──query
    │   ├──cs_dn_n1x476251.gctx
    │   ├──cs_n1x476251.gctx
    │   ├──cs_up_n1x476251.gctx
    │   ├──dn.gmt
    │   ├──up.gmt
    │   ├──leadf_dn_n1x476251.gctx
    │   └──leadf_up_n1x476251.gctx
    └──gutc
        ├──cs_sig.gctx
        ├──ns_pcl_cell.gctx
        ├──ns_pcl_summary.gctx
        ├──ns_pert_cell.gctx
        ├──ns_pert_summary.gctx
        ├──ns_sig.gctx
        ├──ps_pcl_cell.gctx
        ├──ps_pcl_summary.gctx
        ├──ps_pert_cell.gctx
        ├──ps_pert_summary.gctx
        └──query_info.txt
```
