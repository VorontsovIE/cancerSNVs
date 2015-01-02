## Tools and processing workflow.

### bin/compare_actual_site_motifs/:
Extract site sequences(FASTA), prepare PCMs from extracted FASTA sequences, compare them and draw logos for inconsistent motif pairs.

```sh
ruby bin/compare_actual_site_motifs/collect_motif_sites.rb sites.txt ./fasta/
ruby bin/compare_actual_site_motifs/collect_motif_sites.rb sites2.txt ./fasta2/
ls fasta/* | ruby bin/compare_actual_site_motifs/fasta2pcm.rb --stdin-filelist
ls fasta2/* | ruby bin/compare_actual_site_motifs/fasta2pcm.rb --stdin-filelist
ls fasta/pcm/* | sequence_logo --logo-folder fasta/logo/
ls fasta2/pcm/* | sequence_logo --logo-folder fasta2/logo/
ruby bin/compare_actual_site_motifs/compare_motif_pairs.rb fasta/pcm/ fasta2/pcm/
```

Here `sites.txt` and `sites2.txt` are two site lists to compare. E.g. sites before vs after deleting repeats. Or sites in cancer vs sites in shuffle.
