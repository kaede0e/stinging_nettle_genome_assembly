conda activate BUSCO_env

# Run BUSCO job:
busco --offline \
--in /~/your_protein.faa \
--out Nettle_Round8_H1_protein_BUSCO_output_with_augustus \
--lineage_dataset /~/busco_downloads/lineages/eudicots_odb10 \
--mode proteins \
--augustus \
--cpu 20

