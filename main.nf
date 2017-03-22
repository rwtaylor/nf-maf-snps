#!/usr/bin/env nextflow

sample_groups = []
samples = []
file("${params.sample_groups}").readLines().each{line ->
  (sampleID, populationID) = line.split(/\t/)
  sample_groups.push([ "p":populationID, "s":sampleID ])
  samples.push([sampleID])
}

sample_groups = sample_groups.groupBy({it -> it.p}).collectEntries{[it.key, it.value.s]}
sample_groups = sample_groups.collect { key, value -> [key, value] }
sample_groups.each{a -> 
  a[1].removeAll(params.excludesamples)
}
sample_groups.removeAll({ it[1].empty })
sample_groups_file = file(params.sample_groups)

println(sample_groups)
println(samples)

input_vcf = Channel.fromPath(params.vcf_file).map { file -> [file.baseName, file] }
input_vcf = input_vcf.view()
input_vcf.into{input_vcf_plink; input_vcf_csvcf}

process Plink_bed {
  publishDir 'outputs/plink', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(vcf) from input_vcf_plink

  output:
  set prefix, file("*.bed"), file("*.bim"), file("*.fam") into plink_bed

  script:
  remove_file = file("temp-remove.txt")
  remove_file.text = "temp temp\n"
  params.excludesamples.each{
    remove_file.append("${it} ${it}\n")
  }
  """
  /usr/local/bin/plink --make-bed --remove $remove_file --vcf $vcf --allow-extra-chr --out temp
  /usr/local/bin/plink --make-bed --set-missing-var-ids @:#\\\$1,\\\$2 --allow-extra-chr --bfile temp --out $prefix
  rm temp.bed temp.bim temp.fam
  """
}

plink_bed.into{plink_bed_for_pruning; plink_bed_ldak; plink_bed}

process Plink_ld_pruning {
  publishDir 'outputs/plink', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(bed), file(bim), file(fam) from plink_bed_for_pruning

  output:
  set val("${prefix}-ldp"), file("*.bed"), file("*.bim"), file("*.fam") into plink_pruned_bed
  set val("${prefix}-ldp"), file("*.prune.in"), file("*.prune.out") into plink_pruned
  
  """
  /usr/local/bin/plink --indep 50 5 2 --allow-extra-chr --bed $bed --bim $bim --fam $fam --out $prefix
  /usr/local/bin/plink --make-bed --extract ${prefix}.prune.in --allow-extra-chr --bed $bed --bim $bim --fam $fam --out ${prefix}.ldpruned
  """
}

// Combine pruned and non-pruned genotypes
plink_bed = plink_bed.mix(plink_pruned_bed)
plink_bed.into{plink_bed_flat; plink_bed_stats}

process Plink_traw {
  publishDir 'outputs/plink', mode: 'copy'
  tag {prefix}
  cpus 1
  memory 4.GB
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(bed), file(bim), file(fam) from plink_bed_flat

  output:
  set prefix, file("*.traw") into plink_traw

  """
  /usr/local/bin/plink --recode A-transpose --allow-extra-chr --bed ${bed} --bim ${bim} --fam ${fam} --out ${prefix}
  """
}

plink_traw.into { plink_traw_allpoly; plink_traw_diff }

process AllPoly_bed {
  publishDir 'outputs/all_poly-bed', mode: 'copy'
  tag {prefix}
  cpus 16
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 16.GB}
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(traw) from plink_traw_allpoly

  output:
  file("*.bed") into allpoly_beds

  """
  all_poly_snps.R $task.cpus jacksoni,sumatrae,altaica,tigris $sample_groups_file $traw $prefix
  """
}

process Differentiating_bed {
  publishDir 'outputs/differentiating-bed', mode: 'copy'
  tag {prefix}
  cpus 16
  memory { task.exitStatus == 137 ? (task.attempt > 2 ? 64.GB: 32.GB) : 16.GB}
  time 1.h
  errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
  maxRetries 7
  maxErrors '-1'

  input:
  set prefix, file(traw) from plink_traw_diff

  output:
  file("*.bed") into diff_beds

  """
  differentiating_snps.R $task.cpus jacksoni,sumatrae,altaica,tigris $sample_groups_file $traw $prefix
  """
}

workflow.onComplete {
  println "Pipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
