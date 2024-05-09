version 1.0

# example WDL style guidelines:
# https://broadinstitute.github.io/warp/docs/contribution/contribute_to_warp/wdl_task_runtime_style/

# optimized for use in Terra with GCP C3 machines described here:
# https://cloud.google.com/compute/all-pricing (as of 2023/10/17)


workflow TotalRecall {

  input {
    String sample_name
    File case_bam
    File case_bam_index
    File control_bam
    File control_bam_index
    String docker_image
    Int preemptible_tries = 3
    Boolean debug = false
  }

  call TotalRecall_pp_bam as pp_case {
    input:
      bam = case_bam,
      bam_index = case_bam_index,
      sample = "case",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_sort_clipped_bam as sort_bam_case {
    input:
      bam = pp_case.clip_bam,
      sample = "case",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_clip_tsv as case_clip {
    input:
      clip_bam = sort_bam_case.sorted_bam,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_drp_names {
    input:
      case_disc_bam = pp_case.disc_bam,
      final_input_p1 = TotalRecall_integration_sites.final_input_p1,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_drp_fq {
    input:
      case_disc_bam = pp_case.disc_bam,
      case_drp_names_zst = TotalRecall_drp_names.case_drp_names_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_drp_aln {
    input:
      case_disc_fq_zst = TotalRecall_drp_fq.case_disc_fq_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_drp_case {
    input:
      case_drp_subset_bam = TotalRecall_drp_fq.disc_subset_bam,
      ri_tsv_zst = TotalRecall_drp_aln.case_disc_ri_zst,
      case_disc_fq_zst = TotalRecall_drp_fq.case_disc_fq_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_pp_bam as pp_ctrl {
    input:
      bam = control_bam,
      bam_index = control_bam_index,
      sample = "control",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_sort_clipped_bam as sort_bam_ctrl {
    input:
      bam = pp_ctrl.clip_bam,
      sample = "control",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_clip_tsv as control_clip {
    input:
      clip_bam = sort_bam_ctrl.sorted_bam,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_case_1p {
    input:
      clip_tsv_zst = case_clip.clip_tsv_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_aln_1p {
    input:
      case_1p_fq_zst = TotalRecall_case_1p.case_1p_fq_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_breakpoints_2p {
    input:
      case_clip_tsv_zst = case_clip.clip_tsv_zst,
      case_1p_other_tsv_zst = TotalRecall_aln_1p.case_1p_other_tsv_zst,
      case_polya_zst = TotalRecall_case_1p.case_polya_zst,
      case_mask_tsv = TotalRecall_case_1p.case_mask_tsv,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_aln_2p {
    input:
      case_clipped_fa_zst = TotalRecall_breakpoints_2p.case_clipped_fa_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_integration_sites {
    input:
      case_bp_byid_zst = TotalRecall_breakpoints_2p.case_bp_byid_zst,
      case_clip_annot_tsv_zst = TotalRecall_aln_2p.case_clip_annot_tsv_zst,
      case_md = pp_case.md,
      case_clip_tsv_zst = case_clip.clip_tsv_zst,
      control_clip_tsv_zst = control_clip.clip_tsv_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_clipped_clusters as cc_case{
    input:
      case_ins_nbp_zst = TotalRecall_integration_sites.case_ins_nbp_zst,
      clip_tsv_zst = case_clip.clip_tsv_zst,
      case_mask_tsv = TotalRecall_case_1p.case_mask_tsv,
      sample = "case",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_clipped_clusters as cc_ctrl{
    input:
      case_ins_nbp_zst = TotalRecall_integration_sites.case_ins_nbp_zst,
      clip_tsv_zst = control_clip.clip_tsv_zst,
      case_mask_tsv = TotalRecall_case_1p.case_mask_tsv,
      sample = "control",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_depth_info as doc_case{
    input:
      bam = case_bam,
      bam_index = case_bam_index,
      case_md = pp_case.md,
      final_input_p1 = TotalRecall_integration_sites.final_input_p1,
      sample = "case",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_depth_info as doc_ctrl{
    input:
      bam = control_bam,
      bam_index = control_bam_index,
      case_md = pp_case.md,
      final_input_p1 = TotalRecall_integration_sites.final_input_p1,
      sample = "control",

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  call TotalRecall_final_step {
    input:
      case_bam = case_bam,
      case_bam_index = case_bam_index,
      control_bam = control_bam,
      control_bam_index = control_bam_index,
      case_md = pp_case.md,
      control_md = pp_ctrl.md,
      final_input_p1 = TotalRecall_integration_sites.final_input_p1,
      cc_case_tar = cc_case.cc_tar_zst,
      cc_ctrl_tar = cc_ctrl.cc_tar_zst,
      doc_case_tar = doc_case.doc_tar_zst,
      doc_ctrl_tar = doc_ctrl.doc_tar_zst,
      case_disc_tsv_zst = TotalRecall_drp_case.case_disc_tsv_zst,

      docker=docker_image,
      preemptible_tries=preemptible_tries,
      debug = debug
  }

  output {
    File results_tar = TotalRecall_final_step.results_tar
  }
}

task TotalRecall_pp_bam {
  # approximate run time: 4h

  input {
    File bam
    File bam_index
    String sample

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 3
    Int memory_gb = 3
  }
  Float bam_size = size(bam, "GiB")
  Int disk_size_gb = ceil((1.6 * bam_size) + 20)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{bam} ~{sample}.bam
    mv ~{bam_index} ~{sample}.bam.bai
    touch ~{sample}.bam.bai

    snakemake --cores $(nproc)  -s /opt/totalrecall/Snakefile -p --nt \
      --set-threads split_case_bam=$(nproc) split_ctrl_bam=$(nproc) -- ~{sample}.clip.bam
    snakemake --cores $(nproc) -s /opt/totalrecall/Snakefile -p --nt ~{sample}.genome
    snakemake --cores $(nproc) -s /opt/totalrecall/Snakefile -p --nt ~{sample}.header.sam

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      df -h . || true
    fi

    tar -c ~{sample}.nbases ~{sample}.rlen ~{sample}.tlen ~{sample}.genome ~{sample}.header.sam | zstd -q > md.tar.zst

    mv ~{sample}.clip.bam clip.bam
    [ ! -f "~{sample}.disc.bam" ] && touch "~{sample}.disc.bam"  # create dummy output for control
    mv ~{sample}.disc.bam disc.bam

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} SSD"
  }

  output {
    File md = "md.tar.zst"
    File clip_bam = "clip.bam"
    File disc_bam = "disc.bam"
  }
}  # end task TotalRecall_pp_bam

task TotalRecall_sort_clipped_bam {
  # approximate run time: 1h

  input {
    File bam
    String sample

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 4
    Int memory_gb = 15
  }
  Float bam_size = size(bam, "GiB")
  Int disk_size_gb = ceil((6 * bam_size) + 40)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{bam} ~{sample}.clip.bam

    if [ ~{memory_gb} -lt 4 ]
    then
      echo "need at least 4GiB of RAM" 1>&2
      exit 1
    fi

    # set aside 3GiB and use the rest of memory for sorting BAM
    sort_bam_mem_mb=$(printf "%d %d\n" ~{memory_gb} $(nproc) | awk '{print int(($1 - 3) * 1024 / $2)}')

    snakemake --cores $(nproc)  -s /opt/totalrecall/Snakefile -p --nt \
      --set-threads sort_clipped_bam=$(nproc) --config bam_sort_mb=$sort_bam_mem_mb -- \
      ~{sample}.clip.namesorted.bam

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      df -h . || true
    fi

    mv ~{sample}.clip.namesorted.bam clip.namesorted.bam

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} SSD"
  }

  output {
    File sorted_bam = "clip.namesorted.bam"
  }
}  # end task TotalRecall_sort_clipped_bam

task TotalRecall_clip_tsv {
  # approximate run time: 4h

  input {
    File clip_bam

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 1
    Int memory_gb = 4
  }
  Float bam_size = size(clip_bam, "GiB")
  Int disk_size_gb = ceil((5 * bam_size) + 50)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{clip_bam} this.clip.namesorted.bam

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/breakpoints.snakefile -p this.clip.tsv

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm this.clip.tsv

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File clip_tsv_zst = "this.clip.tsv.zst"
  }
}  # end task TotalRecall__clip_tsv

task TotalRecall_drp_names {
  # approximate run time:

  input {
    File case_disc_bam
    File final_input_p1

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 1
    Int memory_gb = 6
  }
  Float bam_size = size(case_disc_bam, "GiB")
  Int disk_size_gb = ceil((10 * bam_size) + 20)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    ln -s ~{case_disc_bam} case.disc.bam
    tar -xf ~{final_input_p1}

    snakemake --cores $(nproc) --set-threads drp_names=$(nproc) \
      -s /opt/totalrecall/snakemake/drp.snakefile -p case.drp_names

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T$(nproc) case.drp_names

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # worst case scenario: fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File case_drp_names_zst = "case.drp_names.zst"
  }
}  # end task TotalRecall_drp_names

task TotalRecall_drp_fq {
  # approximate run time: 3h

  input {
    File case_disc_bam
    File case_drp_names_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 4
  }
  Float bam_size = size(case_disc_bam, "GiB") + 4 * size(case_drp_names_zst, "GiB")
  Int disk_size_gb = ceil((10 * bam_size) + 20)
  Int memory_gb = ceil(10 * size(case_drp_names_zst, "GiB")) + 10

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    ln -s ~{case_disc_bam} case.disc.bam
    mv ~{case_drp_names_zst} case.drp_names.zst
    unzstd -q --rm case.drp_names.zst

    if [ ~{memory_gb} -lt 6 ]
    then
      echo "need at least 6GiB of RAM" 1>&2
      exit 1
    fi

    # set aside 4GiB and use the rest of memory for sorting BAM
    sort_bam_mem_mb=$(printf "%d %d\n" ~{memory_gb} $(nproc) | awk '{print int(($1 - 4) * 1024 / $2)}')

    snakemake --cores $(nproc) --set-threads discordant_reads_fastq=$(nproc) drp_subset_bam=$(nproc) \
      --config bam_sort_mb=$sort_bam_mem_mb --nt -s /opt/totalrecall/snakemake/drp.snakefile \
      -p case.disc.fq

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T$(nproc) case.disc.fq

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # worst case scenario: fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} SSD"
  }

  output {
    File case_disc_fq_zst = "case.disc.fq.zst"
    File disc_subset_bam = "case.drp.subset.bam"
  }
}  # end task TotalRecall_drp_fq

task TotalRecall_drp_aln {
  # approximate run time: 3h

  input {
    File case_disc_fq_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 20
    Int memory_gb = 10
  }
  Float fq_size = size(case_disc_fq_zst, "GiB")
  Int disk_size_gb = ceil((25 * fq_size) + 20)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_disc_fq_zst} case.disc.fq.zst
    unzstd -q --rm case.disc.fq.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/drp.snakefile -p case.repeatinfo.tsv

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T$(nproc) case.repeatinfo.tsv

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # worst case scenario: fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} SSD"
  }

  output {
    File case_disc_ri_zst = "case.repeatinfo.tsv.zst"
  }
}  # end task TotalRecall_drp_aln

task TotalRecall_drp_case {
  # approximate run time: 3h

  input {
    File case_drp_subset_bam
    File ri_tsv_zst
    File case_disc_fq_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 4
    Int memory_gb = 20
  }
  Float input_size = size(case_drp_subset_bam, "GiB") + size(ri_tsv_zst, "GiB")
  Int disk_size_gb = ceil((8 * input_size) + 20)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    ln -s ~{case_drp_subset_bam} case.drp.subset.bam
    mv ~{ri_tsv_zst} case.repeatinfo.tsv.zst
    mv ~{case_disc_fq_zst} case.disc.fq.zst
    unzstd -q --rm case.repeatinfo.tsv.zst
    unzstd -q --rm case.disc.fq.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p \
      --set-threads merge_discordant_repeat_info=$(nproc) -- case.disc.tsv

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T$(nproc) case.disc.tsv

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # worst case scenario: fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File case_disc_tsv_zst = "case.disc.tsv.zst"
  }
}  # end task TotalRecall_drp_case

task TotalRecall_case_1p {
  # approximate run time: 3h

  input {
    File clip_tsv_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int cpu = 1
    Int memory_gb = 8
  }
  Float zst_size = size(clip_tsv_zst, "GiB")
  Int disk_size_gb = ceil((10 * zst_size) + 10)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{clip_tsv_zst} case.clip.tsv.zst
    unzstd -q --rm case.clip.tsv.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/breakpoints.snakefile -p --nt case.1p.fq

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T4 case.clip.tsv
    zstd -q --rm -T4 case.1p.polya.tsv
    zstd -q --rm -T4 case.1p.fq

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File case_polya_zst = "case.1p.polya.tsv.zst"
    File case_1p_fq_zst = "case.1p.fq.zst"
    File case_mask_tsv = "case.mask.tsv"
  }
}  # end task TotalRecall_case_1p_tsv

task TotalRecall_aln_1p {
  # approximate run time: 3h

  input {
    File case_1p_fq_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 10
    Int cpu = 30
  }
  Float fq_size = size(case_1p_fq_zst, "GiB")
  Int disk_size_gb = ceil((25 * fq_size) + 10)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_1p_fq_zst} case.1p.fq.zst
    unzstd -q --rm case.1p.fq.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/lastal.snakefile -p case.1p.other.tsv

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T4 case.1p.other.tsv

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # in the worst case scenario we will fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    # LAST alignment done here is well parallelizable, does not need RAM but needs CPU
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File case_1p_other_tsv_zst = "case.1p.other.tsv.zst"
  }
}  # end task TotalRecall_aln_1p

task TotalRecall_breakpoints_2p {
  # approximate run time: 3h

  input {
    File case_clip_tsv_zst
    File case_1p_other_tsv_zst
    File case_polya_zst
    File case_mask_tsv

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 8
    Int cpu = 1
  }
  Float tsv_zst_size = size(case_clip_tsv_zst, "GiB")
  Int disk_size_gb = ceil((10 * tsv_zst_size) + 10)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_clip_tsv_zst} case.clip.tsv.zst
    mv ~{case_1p_other_tsv_zst} case.1p.other.tsv.zst
    mv ~{case_polya_zst} case.1p.polya.tsv.zst
    mv ~{case_mask_tsv} case.mask.tsv
    unzstd -q --rm case.clip.tsv.zst
    unzstd -q --rm case.1p.other.tsv.zst
    unzstd -q --rm case.1p.polya.tsv.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/breakpoints2p.snakefile --nt -p case.clipped.fa

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm case.breakpoints.byid.tsv
    zstd -q --rm case.clipped.fa

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # in the worst case scenario we will fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File case_bp_byid_zst = "case.breakpoints.byid.tsv.zst"
    File case_clipped_fa_zst = "case.clipped.fa.zst"
  }
}  # end task TotalRecall_breakpoints_2p

task TotalRecall_aln_2p {
  # approximate run time: 3h

  input {
    File case_clipped_fa_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 10
    Int cpu = 30
  }
  Float fa_zst_size = size(case_clipped_fa_zst, "GiB")
  Int disk_size_gb = ceil((30 * fa_zst_size) + 50)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_clipped_fa_zst} case.clipped.fa.zst
    unzstd -q --rm case.clipped.fa.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/lastal.snakefile -p case.clip_annot.tsv

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    zstd -q --rm -T20 case.clip_annot.tsv

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    # in the worst case scenario we will fit into a high CPU/low mem C3 machine with 22CPUs/44G RAM
    # LAST alignment done here is well parallelizable, does not need RAM but needs CPU
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File case_clip_annot_tsv_zst = "case.clip_annot.tsv.zst"
  }
}  # end task TotalRecall_aln_2p

task TotalRecall_integration_sites {
  # approximate run time: 2h

  input {
    File case_bp_byid_zst
    File case_clip_annot_tsv_zst
    File case_md

    File case_clip_tsv_zst
    File control_clip_tsv_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 8
    Int cpu = 1
  }
  Float input_size = size(case_bp_byid_zst, "GiB") + size(case_clip_tsv_zst, "GiB") + size(control_clip_tsv_zst, "GiB")
  Int disk_size_gb = ceil((7 * input_size) + 50)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_bp_byid_zst} case.breakpoints.byid.tsv.zst
    mv ~{case_clip_annot_tsv_zst} case.clip_annot.tsv.zst
    mv ~{case_md} case_md.tar.zst
    mv ~{case_clip_tsv_zst} case.clip.tsv.zst
    mv ~{control_clip_tsv_zst} control.clip.tsv.zst
    unzstd -q --rm case.breakpoints.byid.tsv.zst
    unzstd -q --rm case.clip_annot.tsv.zst
    unzstd -q --rm case.clip.tsv.zst
    unzstd -q --rm control.clip.tsv.zst
    tar -xf case_md.tar.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt is.bed
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt case.ins.namesorted.tsv
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt aux/case.low_complexity
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt aux/control.clipped.read

    tar -c aux/*  case.ins.cand.tsv case.ins.namesorted.tsv is.bed drp.bed | zstd -q > final_input_p1.tar.zst
    zstd -q --rm case.ins.nbp.tsv

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} SSD"
  }

  output {
    File final_input_p1 = "final_input_p1.tar.zst"
    File case_ins_nbp_zst = "case.ins.nbp.tsv.zst"
  }
}  # end task TotalRecall_integration_sites

task TotalRecall_clipped_clusters {
  # approximate run time: 90 min
  input {
    File case_ins_nbp_zst
    File clip_tsv_zst
    File case_mask_tsv
    String sample

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 8
    Int cpu = 1
  }
  Float input_size = size(clip_tsv_zst, "GiB") + size(case_ins_nbp_zst, "GiB")
  Int disk_size_gb = ceil((10 * input_size) + 50)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_ins_nbp_zst} case.ins.nbp.tsv.zst
    mv ~{clip_tsv_zst} ~{sample}.clip.tsv.zst
    mv ~{case_mask_tsv} case.mask.tsv
    unzstd -q --rm case.ins.nbp.tsv.zst
    unzstd -q --rm ~{sample}.clip.tsv.zst
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt aux/~{sample}.read.cluster

    tar -c aux/* | zstd -q > clipped_clusters.tar.zst

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File cc_tar_zst = "clipped_clusters.tar.zst"
  }
}  # end task TotalRecall_clipped_clusters

task TotalRecall_depth_info {
  # approximate run time: 3h

  input {
    File bam
    File bam_index
    File case_md
    File final_input_p1
    String sample

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 20
    Int cpu = 1
  }
  Float input_size = size(bam, "GiB")
  Int disk_size_gb = ceil((1.5 * input_size) + 50)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{bam} ~{sample}.bam
    mv ~{bam_index} ~{sample}.bam.bai
    tar -xf ~{case_md}
    tar -xf ~{final_input_p1}
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/depth.snakefile -p --nt aux/doc_~{sample}.tsv

    tar -c aux/doc*.tsv | zstd -q > doc.tar.zst

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File doc_tar_zst = "doc.tar.zst"
  }
}  # end task TotalRecall_depth_info

task TotalRecall_final_step {
  # approximate run time: 3h

  input {
    File case_bam
    File case_bam_index
    File control_bam
    File control_bam_index
    File case_md
    File control_md
    File final_input_p1
    File cc_case_tar
    File cc_ctrl_tar
    File doc_case_tar
    File doc_ctrl_tar
    File case_disc_tsv_zst

    String docker
    Int preemptible_tries
    Boolean debug
    Int memory_gb = 8
    Int cpu = 1
  }
  Float bam_size = size(case_bam, "GiB") + size(control_bam, "GiB") + 3 * size(case_disc_tsv_zst, "GiB")
  Int disk_size_gb = ceil((1.4 * bam_size) + 50)

  command <<<
    export LC_COLLATE=C

    MON_PID=0
    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      printf "working dir: %s\n" $(pwd) || true
      printf "requested %s GiB disk, %s GiB RAM\n" ~{disk_size_gb} ~{memory_gb} || true
      df -h . || true
      free -h || true
      lscpu || true
      cat /proc/cpuinfo || true
      cromwell_monitoring_script.sh > monitoring.log &
      MON_PID="$!"
    fi

    mv ~{case_bam} case.bam
    mv ~{case_bam_index} case.bam.bai
    mv ~{control_bam} control.bam
    mv ~{control_bam_index} control.bam.bai
    mv ~{case_disc_tsv_zst} case.disc.tsv.zst
    tar -xf ~{case_md}
    tar -xf ~{control_md}
    tar -xf ~{final_input_p1}
    tar -xf ~{cc_case_tar}
    tar -xf ~{cc_ctrl_tar}
    tar -xf ~{doc_case_tar}
    tar -xf ~{doc_ctrl_tar}
    unzstd -q --rm case.disc.tsv.zst

    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt aux/case.pe.L.tsv
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/penultimate.snakefile -p --nt aux/case.pe.R.tsv
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/final.snakefile -p --nt all.vcf.gz
    snakemake --cores $(nproc) -s /opt/totalrecall/snakemake/final.snakefile -p --nt results.tgz

    if [ ~{true="1" false="0" debug} -gt 0 ]
    then
      date || true
      df -h || true
      du -sch * || true
    fi

    if [ $MON_PID -ne 0 ]
    then
      kill $MON_PID && sleep 20
      cat monitoring.log
    fi
  >>>

  runtime {
    docker: docker
    preemptible: preemptible_tries
    memory: "${memory_gb} GiB"
    cpu: cpu
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File results_tar = "results.tgz"
  }
}  # end task TotalRecall_final_step
