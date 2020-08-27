#/bin/sh
# gencat.sh: extract microbial protein records from the RefSeq catalog
# 20190830 updated for RefSeq release 95
# 20190716 updated for RefSeq release 201

# zcat: decompress the release catalog
# awk: save the tax_id, accession, and FTP directory columns
#      (specifically, we drop the taxon name because e.g. Staphylococcus aureus FP_N5208 gives a false match in the regex)
# grep: extract the microbial protein records (148851653 records, 5.5 GB)
zcat RefSeq-release201.catalog.gz | awk -F$'\t' '{print $1,$3,$4}' | grep "P_.*bacteria\|P_.*archaea\|P_.*viral" > microbial.protein.tmp
# save only the accession and taxid columns (3.0 GB)
awk '{print $2,$1}' microbial.protein.tmp > accession.taxid.tmp
rm microbial.protein.tmp
# sort the file on accession so that it can be used with the unix 'join' command
# NOTE: for join to work, -n for numeric sort is *not* used here
sort accession.taxid.tmp > accession.taxid.txt
rm accession.taxid.tmp

# also make a list of unique taxids 20190830
awk '{print $2}' accession.taxid.txt > taxid.tmp
sort -u taxid.tmp > unique_taxids.txt
# also for the WP accessions
awk '{print $2}' release201.MultispeciesAutonomousProtein2taxname > taxid.tmp
sort -u taxid.tmp > WP_unique_taxids.txt
rm taxid.tmp
