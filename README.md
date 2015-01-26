## Workflow
Workflow is following:
* Get cancer SNVs
* Remove and recreate location type markup: promoter and intronic regions
* Filter SNVs to obtain regulatory SNVs only (located either in promoter or in intronic regions). This SNV collection is used for the rest of computations. We don't use non-regulatory SNVs anymore.
* Generate sequences with substitutions for site extraction:
* * Extract sequences around cancer SNVs (regulatory, you remember) from genome.
* * Generate several copies (e.g. 10; more helps in further fitting, but linearly slows down site extraction, which is the slowest step in workflow) of shuffled sequences. Shuffling preserves 1-bp context around SNV, so distribution of mutation contexts is the same as in cancer.
* * Take random sequences from genome (not equal to cancer positions) and put random mutations into them, so that mutation types distribution is the same as in cancer. Better to have several times more mutations than in cancer (e.g. 10 times more)
* (parallelizable task, use several cores on several servers because it's a bottleneck. It can take about 1.5 days using 8 CPU cores for multiplication factor 10 as we used). Run PerfectosAPE in order to obtain sites in all sequences collected in a previous step (cancer, genomic random and shuffling random).
* If necessary, take subset of sites for a certain mutation type, such as TpCpN or NpCpG (on either strand). Because random sequences preserve mutation context distribution, this step is safe. If one need the only mutation context, now and further (not our case), he can make this filtering before sequence extraction and generation.
* Make fitting of random sites to real cancer sites. Sites around somatic mutations and around randomly created mutations can present different P-value distributions, we sample sites of random mutations so that each site's P-value distribution for each different mutation context was proportional to a corresponding in-cancer distribution.
* Calculate statistics of sites total and of disrupted sites in cancer substitutions and in choosen random substitutions. 
* Evaluate Fischer exact test for disrupted sites (across all sites) in cancer vs the same value in random mutations.
* Now we can get rate of site disruption in cancer vs rate of site disruption  in random mutations and calculate their ratio. Significance of this difference is given by a P-value from Fisher exact test (we should make multiple comparison correction such as Bonferony/Holms/FDR).
* Filter only motifs whose sites are disrupted significantly more often than in random, choose only A-,B,-C- quality motifs.
* Repeat the same computations for several random mutation sets. Intersect results in order to get reliable answer.

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
