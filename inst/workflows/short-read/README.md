# TEi short-read WDL workflow

This WDL 1.0 workflow packages the original three-stage TEi short-read design:
extract terminal soft clips from a host-aligned BAM, align them to a TE
reference, and report supported insertion breakpoints.

The execution environment must provide R, TEi, Rbowtie2, and a WDL 1.0 runner
such as Cromwell. The host BAM is read sequentially and does not require an
index. The TE FASTA is indexed inside the alignment task.

Copy `inputs.example.json`, replace the two absolute input paths, and run:

```sh
java -jar cromwell.jar run TEi-short-read.wdl -i inputs.json
```

The workflow was developed and is maintained by Yihan Xiao. The underlying TEi
short-read R functions were co-developed by Yihan Xiao and Tao Chen.
