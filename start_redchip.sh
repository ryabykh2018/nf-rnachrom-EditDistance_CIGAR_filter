cd /home/ryabykh2018/storage/EditDistance_CIGAR_filter/

# for file in SRR17331253_K562_CTCF_RedChIP_rep2_Unique_RNA SRR17331254_K562_CTCF_RedChIP_rep2_Unique_RNA SRR17331268_K562_CTCF_RedChIP_rep1_Unique_RNA; do #SRR17331267_K562_CTCF_RedChIP_rep1_Unique_RNA    
#     /usr/bin/time -v python EditDistance_CIGAR_filter.py "NM + N_softClipp_bp" 2 2 0 0 "yes" "explorer" "ATA, not iMARGI" "$file.tab.rc" "/mnt/smb_share/rnachrom/contacts/redchip/Align_bam/" "/home/ryabykh2018/storage/EditDistance_CIGAR_filter/redchip/$file/"
# done

# for file in SRR9201801_mESC_1FA_n2_Unique_RNA SRR9201803_mESC_1FA_n3_Unique_RNA SRR9201805_mESC_2FA_n1_Unique_RNA SRR9201807_mESC_2FA_n2_Unique_RNA SRR9201809_mESC_2FA_n3_Unique_RNA; do #SRR9201799_mESC_1FA_n1_Unique_RNA.tab.rc
#     /usr/bin/time -v python EditDistance_CIGAR_filter.py "NM + N_softClipp_bp" 2 2 0 0 "yes" "explorer" "ATA, not iMARGI" "$file.tab.rc" "/mnt/smb_share/rnachrom/contacts/radicl/Align_bam/" "/home/ryabykh2018/storage/EditDistance_CIGAR_filter/radicl/$file/"
# done
