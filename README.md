<div align="center">
  <center><h1>KhufuPAN</h1></center>
  <img width="500" alt="Image" src="https://github.com/user-attachments/assets/52299e9b-44ab-485b-9e9d-4735f32af7bf" />
</div>

KhufuPAN is an open-source pipeline for highly efficient genotyping using pangenomes.

## Requirements
### KhufuPAN
- **KhufuEnv:** https://github.com/w-korani/KhufuEnv
- **VG:** https://github.com/vgteam/vg
- **KMC:** https://github.com/refresh-bio/KMC


## Implementation & Getting Help

1. Bootstrap of the pangenome graph
   ```
   bootstrap.sh -h
   ```
2. Processing of fastq libraries
   ```
   LibraryProc.sh -h
   ```
3. Combining a subset and filtering for sequencing depth
   ```
   CombineDS.sh -h
   ```
4. Combining subsets and filtering for the pouplation structure
   ```
   KhufuEnvHelp
   ```

## Example code
```
bootstrap.sh -gfa TestProc.gfa  -t 4
```
```
LibraryProc.sh -gfa refPROC/TestProc.gfa -r1 fastqs/sample1.fq.gz -id S01 -t 4
LibraryProc.sh -gfa refPROC/TestProc.gfa -r1 fastqs/sample2.fq.gz -id S02 -t 4
LibraryProc.sh -gfa refPROC/TestProc.gfa -r1 fastqs/sample3.fq.gz -id S03 -t 4
```
```
CombineDS.sh -gam gams -gfa refPROC/TestProc.gfa -min 1 -max 10 -o TestSet1.panmap -t 4
```
```
panmapFilterMissingSample TestSet1.panmap 0.9
panmapFilterMissingVariant TestSet1.Smiss0.9.panmap 0.75
panmapFilterMAF TestSet1.Smiss0.9.miss0.75.panmap 0.05
panmapGetDiAlleles TestSet1.Smiss0.9.miss0.75.MAF0.05.panmap
panmapGetDiAlleles TestSet1.Smiss0.9.miss0.75.MAF0.05.DiAllelic.panmap
panmapGetHomoPolymorphic TestSet1.Smiss0.9.miss0.75.MAF0.05.DiAllelic.DiAllelic.panmap
panmapGetSV TestSet1.Smiss0.9.miss0.75.MAF0.05.DiAllelic.DiAllelic.HomoPolymorphic.panmap
panmapAlleleFreqStats TestSet1.Smiss0.9.miss0.75.MAF0.05.DiAllelic.DiAllelic.HomoPolymorphic.SV.panmap > stat1
hapmapAlleleTypeFreq TestSet1.Smiss0.9.miss0.75.MAF0.05.DiAllelic.DiAllelic.HomoPolymorphic.SV.panmap > stat2
```

<div align="center">
  <img width="500" alt="Image" src="https://github.com/user-attachments/assets/5ef2d8ff-e97b-4b97-a6f6-b63ba1a0b656" />
</div>


## Citation
Long-Read Low-Pass Sequencing for High-Resolution Trait Mapping
Kendall Lee, Walid Korani, Philip C. Bentz, Sameer Pokhrel, Peggy Ozias-Akins, Alex Harkess, Justin Vaughn, Josh Clevenger
bioRxiv 2025.01.09.632261; doi: https://doi.org/10.1101/2025.01.09.632261
