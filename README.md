![IMB-logo](resources/IMB_logo.png)

# GLOE-Pipe

The **G**enome-wide **L**igation of 3'-**O**H **E**nds sequencing (GLOE-Seq) is a sensitive method to capture DNA strand breaks genomewide at nucleotide resolution.

To facilitate the analysis of the raw GLOE-Seq data, we have developed an easy-to-use, modular and versatile computational pipeline, called **GLOE-Pipe**.

Given that the sequences produced by GLOE-Seq are the reverse-complement of the original captured DNA fragment, GLOE-Pipe assigns the 5'-end of a read to the strand opposite to that on which it mapped. 

Two output modes are available:
- in the **direct** mode, each read represents one unit of signal that is positioned exactly at the 3'-terminal nucleotide of the original captured fragment. This mode would normally be used to visualise strand breaks
- the **indirect** mode is essentially identical, but it positions the signal one nucleotide immediately upstream of the 5'-end of the original captured fragment, which corresponds to the nucleotide immediately 3' of the break. The indirect mode would therefore normally be used to visualise the position of modified or damaged bases (under the premise that the enzymatic treatment used to detect this modification or lesion generates a nick 5’ to the affected nucleotide)

A flowchart of GLOE-Pipe, detailing each one of its steps, can be found [here](https://www.draw.io/?lightbox=1&highlight=0000ff&edit=_blank&layers=1&nav=1&title=GLOEseq_pipeline.html#R7V3rd5u4Ev9rfG73g314%2BvExTuKmvb1pN0lPd%2FdLDwbZZoMRBRwn%2BeuvRgjMQ2Bs83KTds8WhMCSZjTzm9Fo1JMv188fXc1Z%2FQ8byOpJgvHck696kiSqQ4X8AyUvQclEloOCpWsarNKu4N58RaxQYKUb00BeoqKPseWbTrJQx7aNdD9Rprku3iarLbCV%2FFVHW6JMwb2uWdnSH6bhr1ipOJzsHtwgc7liPz2WRsGDuaY%2FLl28sdnv2dhGwZO1Fn6G9dFbaQbexork65586WLsB1fr50tkwbCGIxa8N8t5GjXZRbZf5oUb4fNidTv2rafV3fz25n728%2FWmHxLgSbM2bCx60tAiH5wa5hOMrmUubfpg%2BGsDTZ26dBiiW3K1ZP%2FS1%2BZuuoQ0iX4rLKWD4b%2BEY7%2Fy1xa5EskzS5sjaxoN6SW2sEsryTP6h1TxfBc%2FRkQigzhdYNtnHCUOod2at0IG%2ByL9TnS3MC0r9tHrIfyNPho%2BoTSUp0tXM0wytqliHa9NndwKUMXSPI9dR%2BQVok7GacPI9YRcHz3HihitPiK8Rr77Qqqwp5MJIw2bUiEbbWP8OWRlqxhrygor1NicWEaf3jEHuWD8weeVJ%2F32z%2Fsfm8%2Bbq9d%2FHj49zP72nv7pi%2B%2FM0lVmkQUxySxillsUYZLllqFQAbe8%2FvCvHoS%2F%2Fvk8u7t216Lx8PHlU1%2FO0A4ZROayW%2Bz6K7zEtmZd70pjBIBx2dX5grHDqPIv8v0XRkFt42NSFOMJ9Gz6f8HrA5Xd%2FR17cvXMvkxvXkJaZJgklyAe3rg6KpgirNPQ00KyucjSfPMpqX14JGCvfsMmacqO3BMpSe5xioq%2B5i6Rz95KETJqxvG0VcVzIS6hoPsSewlu%2F44%2F271G73bvGReAK8itDvMWZjEUzkwrbE4NfBOQrahe5Qx2kkaQc%2FUByNkEg4QaAB70PUriC1JBVJxnjnq4lHoX07Vm2jCepoMsk4hTUjgt1hzBryaLA810bEMS3UmWeY5mlynbpAs0kMEuWiSas%2FJ9wLoXQBBptjR9ot4GRH2Qm4%2FfCLGwZ9oYbr58ve5%2FI0NCrucWngPdNM9HZDhm4Uh5rJ6HfsG3LIzI1SB8OiD6Ej%2B9ZDubXzEcX40z5pnukTLeMOQX8wY3DQMqp%2BWxX%2BnT2Wf62py8HWNOYvughWmbvomBaT%2FM4ckfB3Nrp%2FttYc3wqHDxHGKFgXghFpBDmVFYY2NjAeudf%2F%2BQ5xEJbmog%2FzRi2YWd60O3zQXRBVTKa0SIwrw7osMV8OAz0jc%2BjHeMB3Nako%2Bhq4CeYhKLjLPIUx5xkKco1WWnKBytlBoCgr8duCQC1iFIHTofgI9QnY%2FDgsjqV7jDVagWS4%2BhmBkyiWPZhWUngse%2BpCYINkp%2BAC8WHjoVNfLtR7VN1Bhe%2F51AkPtQYwIz7iBkDmqsEhKG02o%2FJmwIEn73kPt1%2Fi%2B434ioBOOb5wLIirQ5fgaJZtrLQKbNsWsgt0%2BKyQ0dopDxqcQDg1SNnsT9a%2FCYGOJIgb9RDUczjOjbUpGC4MGPhrFYQExv4D9znCLJh0WYi9eT%2Fa6UoDRELllC5VMhrZXQOt56QSc8QpqkUzLEPihkakXEitWS8inmz7HxkinMdNQ3wpJP3wqHY1cx%2B2pYYsa%2FJcwRqP%2BFaSGbqPvYp8xSn4fCbHMLOxD%2FnYo64Wlrh9J877fr6tXtt%2B8PFXbKdjZ%2BLnEa71jFFGOda51oIHadY3vBXmaNF%2FAi1qMqGg6FPOFAZSwVcHHwWwxLpJi%2BZDgmBlEAxZm6Zl0w37YPuGMaerottICfwaTWwqLuogX1EE2rBdlK0r0bOf4TiJEDGY8B2eQ2puoPwHfRQtc7viuD70r7%2FDrl8hPzfX7vAHA%2FAKTeEvIpe2Eusw1NPX2HgEdAwIwqon456O4Wu4%2Fg4gF%2FjmZrSyAT%2BLLoTPWyg3gPZjBztLBVwLCVFtbBOxTrmgkazsbwho0QkWiDIkV3NEwo%2BeLByv69YaejjI6DhmiNeB9oGNblmJPG7wChPEAI3ZadAQj7HEAOV6WsSfNNO9AagrNT%2BUF5n84MeDaMPSMM7vfZbLmgesMOPO15YSRhwbfdcgzPQT8vJRycdNkqU%2BuorgbTHh6O%2BShlSvFNC1glWruJIZVoEWTwpLle7qJhXrUsdmHwrWDgc9YHM%2FSoZ%2FRTsMFx8ZNpACcll4M068UzvRwGK1T0tDF0pG%2FJi4nxSowNpwmF60%2FChwB8UNjxR9nW8IZ3r60sFdvKhyi47cr00T2BjPDu1tWcepTeSC1rKasNWso8I6rRJarwBzu7RqUkyKgkP1DJGlX%2Fs2ivXl%2F79y%2F27FJSlrffpde%2BqLQCUSqGGykocGRMUmkHxfBEvEFfvXBd7SVWwQGO8PJZZBSyWhgDNxzHGeLQ%2BuQiaEFJdtoHhiKoMSNK9s9LKuHL4JBfGyJEfRgInVjcLgaOpd5TV9vSMaUhCtAPgnMee2xBvx4IsCAPf%2BmD4FGk2vcqiYBxcpUET7ClFAfzJHyhd1fyMEdd0GBe9mGhcvNJlSZJhlGUjAjkRmAfYT%2BV0STcSMxoNvPlFbLmeNu6NZUVU5HA6wsDQUhGZg6G6niP3KN335BrBghoTxj1fkmn5IjN%2FaJvUoXoy8gqVU6iGEVOBWnkxPceKkOV1B4DWU5tIUm3SzywvpCqP9lTPyWj5fSWlsz302hvdIpMz0y5%2BWhjfMS3fzlo3f%2FzbuJP9I%2FDvng2wc%2BVRigrJaeEwp9MDaxXcAGd3OqilDAYHRWrPjooWD2xiaVeTFiaC8b1EL007CJ6zzeRVBp33SG2MOA4sDrREMJirTwcYikNQax6Y0Xl0bhVSMUfW%2BYIiBnnMbKmZAl4KmODyxyVSWgKg5deEFibhkFFTYoOu5I71nu5l7PlzMW%2BFqxCBFKAKXvSLHVK%2FiMjdgnKQiXNvCT34u6e%2FAfVXf8S2%2BTTmkkJigijbpHn97K75CogvSKIg2TcqZSF01zPjKwcTPuShB5mCF3ew0gDs5mfUTjCdf114wdxPrBXGVjmIiPCyMNfepH0as9XHcgrME2jbphUtpHuuIl2l3E5nr5wFmfYcJvnTFubFvDaDbKeEHy1EQfkSEmZjXJWxjW77BY2oDN8LmRY3cJLr5ucTlrWD3lbx5ZFlEkQdUj5HGn6qhe657vG4DUw8yjcXN4EM%2FMdtMM28Xyba8iVAHiF49Tlj7PUptUmvdmtIs1SWTqVyke5nYbjpFRRhGI3T5Rt5dj6opDiyHpc%2FXNk6CC6abKZknan7iLNp2He2F0TNfBKtcvcXG4hBFAg6Fx%2FBN3zAYC6bYQrwNydpQ3ZrLteHmy3Kk0tDdS7FCCHDNaS3cqX2KM2ZGa9ckzNyjFuxXGrPkapnUXj30ZbcajMH2f5RCpXo63G%2B7TP%2BMT61S5icIdSaRVgHb6IEUvEJNA%2FnRE9yqlLgCd7upEhzbfHAA7%2FzGAH6ebhiEP9LRCHInQQcSjSues9uZTea174DMsKn1GbuEfpkC9GOo4BRmfOAKrQJgOorZgcfAYYHccAUjno2wkO4E%2FCtlfa15rulV9nnxP4QZGFrllWgyvttJWHo4d4aOv5oIf0OvtQ6SB6UM%2FQah4IgpKynIf7BAg3WLEtqTIqK1VqAhblrZpQTmi2DYEOkFavrIyJvRImcHsyPQisfg1LaUh19Bs0z%2BQS2TRbryQsiIW0cVFjcdZBO37umn24nIqT63zkVNrK6aSckjuEcusEuXGRk4gorFL%2BjMsaNq24%2BtLxzepYiXPU%2Fnhodl5Ara67c9xFfVjy5XpV4KSkCpROVYGVeJtVYVLIgqIwVIteqG2x010Yh%2FscXeRYpk51XH%2BB3UdSZJgu5JIFZgz2PYW%2Bx2Z0L%2BnG4do27vI9X20bLjR0Stuq569tO%2BpRUIWy2L9tj0KEyw1iMdF4u9KChkw5zTU9ivG3pg%2FBeejZCYP2wg%2B3toYRAP2wWwcLnoS373wFjzrsoOCRzz7krNxiRl2wqvTOoFOR%2FWnbwSZtUvm0zXv5eqMVOp666%2F%2B0pad26Hhmy4ilaanWQ8vyGl9bS3NkXB2g6G3QP6DTtTV0EffCfMCBw880IGfPAnoh%2F%2BfrDQypTQ93%2BBAYHD0atc8%2FoaMhMMD6fLj58Xts%2FlNaRgEFgYSd3isT2sdd3C0zN5c%2FkgFEwT6ZlPu961vDKuD29M4ZRc1ye207Z7h5C84e3zaZfHFcVnV2Kztz2OyjxFft2Rnv0AIRUujgCfwIq22onLJvX64FyfumX388fLqWTDIlninUWCY6kc71Fz4VII8S5G4ken8TdL%2Fp3Ia09SbWe4n1VJga%2BiPsVgkS9XM9I51Jwwi7MAGY4cB%2FHK3b7stAeVQWxI4lNUwipyFnF6Y4mWR1SX1nbI2zSRIyA1dvBsNI1L3dDIb8s1lbjck77mzWw%2FR8w2ezimLZrQDdwgJRu1uxZcJjbArTW7ja1tBANXbRlqHZD3u7zBYfXk3HSSbY%2FV2MlVTSXM5R81xTJX3qdHUZKyadt8IhURRgwC5ybuASY3z7Bm1tWWgw5cqN8Hmxuh371tPqbn57cz%2F7%2BXrTz%2FciZdjktKNgvyBiZcTZMOcM2nkuqVkGq8yW3lS2KaBwqiiT9Co3yRUPFichRwVcIArpbJNSFiUrBYivci7IZvlOHpZcCjHvd93WMXY83ywPLSt1jd0pjoxhWgPEzukJh3NXSAeWZ4hqa%2BBTe%2B7BP7GDvnfUaV3SN2DS1sFdolRyaiqHJ%2BUvx15Nmq%2BF%2FP12rVfusOSn5Ovs6XLsJKhZ4J8qmI9v9SS29zPFzqJhvbM%2B7CyrY8YldUwV%2BUW5smzUohFrBu6XGU4Zs51cW3kDJmqGPSU5q9D5EGhYE3vmZwXdnUXtk%2FHU47EGu8UUStCNRw%2BgYtzGjqB6QyamNOYksG7Uxsx6yq6wH0Q7a%2BAkD5bJDNh0gdhlLPyJ0gk5XtraeUMkjMKS9pJwcjoJuYEZZ7gPutdEYAb%2F9I1et9ZizuBQzKlJjLsF5HnyqTAXfIytwgCybsGDaa9Dx2LC2A2iGITcYzHzqoWDmTnSu4MHYnZv3ImY9%2FYPe7ZW7qjnjHFByMguHoWLitIfiJ0XCm0LVJXHwleCSOV4QEtw0hv6tSFKGp7pG8%2FHa86UPeloz3zBmhtU3IUgmOQipSqWXOQRBflgqV4mfpg7ihIHVNflVyzSj2%2FYr8gfllaPgTrtsLzjwVLZbAvhCRGto6UIibCjFuiyyxc4HKL8xhCHbiLXvEcqW5eB6O3FjiphhaFHBLpQ126PtemBp1WzEd7APesWXQb0wo0fEnfrRz5Bz27nhyIpqbN%2FeOdINLn3gy8jWt16LgwOMsXITTqnVNXpNo4yzzixckWiqflYucJml9LbG9d6mcL%2BFujgvumz44%2FgLkxlddWfHDj2B8y3iZqab6o8HGQPQhc5IcPyqLYtm2fv52h%2Fm%2Bd0e3G97n%2F9cb9Rrq7%2F%2B3A7%2FLl9Lu8SqWnndWml%2FkDk1RpTN0Q8IJRnuUTvsHMZfXi1wXSR%2Fq6pRfsz8%2BlxdkpaFuWEyJDGcqsqmju0Wanc0dj2A%2BR6dSqQb3PJEylBWFlJWV2BgMkcJr3%2FQ6Nyp1IfYcMVHEnV6ajgYNvrz3hS%2FU4GCP9O23QzIo4zy8qHDo%2BzhkmzB9xJ3efy5JkQ3WTxe9pG%2BPCumaklzqxnOzMv%2BK%2B07pnPtrxg%2FpZq7VlP4kkqBWIHJnE29Pt9ElcyiYXyc7SRyfhWFedkOEzMuUbzWxTka%2B30nIOcqZ2caHezq%2FJT613v5Ogdvpet2VkQOtbbMpFHx%2Bz%2FjuVMLe%2BwLudnaz%2BTeat58lSxTXZoLnn9uTDDqYn26KsHp7UPV%2B%2FLprVP1VeUZlKKs4x4n2rKAgjUpqXJwz8%2BpGJi208J%2BOnwlIDnuTCcTgmohmcCdSkxsHIGuBJ5GwskwGxKE0xfxHI8dRFuxtsnhJuLTjbb4v%2Fn%2FuwsONeHPL4yiXIw5xvWhKnmOhb2eRO%2FfddObqvvQGlQYUfl3sP9PRWEia7kgPBu9XAaiuJLdq48yJX9nehEkx26AeZy5eI19oJ0b42P%2F17G%2F%2B7cIz%2FVrMyvcbanlt3G%2Bpa3hVW6%2FDBkh7ck7EdFzipAsb7M2vn7wrqn7oKzJzup467MBc12CZkJQ7hJQaTne0Wtq6Mt19kjMHZna7QkTd%2FUvFZlnndUrArYEjPBj1tjBKuv%2FkdMGajxfw%3D%3D).

## System Requirements
GLOE-Pipe is designed to run in a high-performance cluster (HPC) environment with a Linux distribution (*e.g.* Debian 9), Bpipe (version 0.9.9.3) as the domain specific language (DSL), and Lmod (version 6.6) as the module system.

## Prerequisites

### Software

| Software | Version used | Link |
| --- | --- | --- |
| `FastQC` | 0.11.5 | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ |
| `Trimmomatic` | 0.36 | http://www.usadellab.org/cms/?page=trimmomatic |
| `Bowtie2` | 2.3.4 | http://bowtie-bio.sourceforge.net/bowtie2/index.shtml |
| `Samtools` | 1.5 | http://samtools.sourceforge.net/ |
| `BEDTools` | 2.25.0 | https://bedtools.readthedocs.io/en/latest/ |
| `bedGraphToBigWig` | 365 | https://github.com/ENCODE-DCC/kentUtils |
| `MACS2` | 2.1.1 | https://github.com/taoliu/MACS |
| `deepTools` | 3.1.0 | https://deeptools.readthedocs.io/en/develop |
| `R` | 3.5.1 | https://www.r-project.org/ |
| `ChIPseeker` | 1.8.0 | https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html |
| `rtracklayer` | 1.42.2 | https://bioconductor.org/packages/release/bioc/html/rtracklayer.html |
| `regioneR` | 1.14.0 | https://bioconductor.org/packages/release/bioc/html/regioneR.html |
| `Cairo` | 1.15-10 | https://cran.r-project.org/web/packages/Cairo/index.html |
| `openxlsx` | 4.1.0 | https://cran.r-project.org/web/packages/openxlsx/index.html |

### Files
- `targets.txt`, which includes sample names
- raw reads (`*.fastq.gz`)

## Before you run GLOE-Pipe

### How to get GLOE-Pipe into your project
This can be done by cloning the most recent development version of GLOE-Pipe from this Git repository

    git clone https://gitlab.com/GPetrosino/GLOE-Pipe.git <project_dir>/GLOEPipe
    cd <project_dir>
    ln -s GLOEPipe/pipelines/GLOEseq/* . 
    ln -s GLOEPipe/modules/GLOEseq/essential.vars.groovy .

### Customise GLOE-Pipe to your needs

Adjust the information found in the following files:

- `gloeseq.pipeline.groovy`, which describes the pipeline steps and the location of the respective modules
- `essential.vars.groovy`, which specifies the main project variables, such as the path to the project directory and the reference genome 
- `targets.txt`, which contains file names and sample comparisons
- `tool.location.groovy` and `bpipe.config`, which specify the paths and resource allocation for each module

## Run GLOE-Pipe

Copy the input FastQ files into `<project_dir>/rawdata folder`.

Using GNU `screen` (for persistence), load the bpipe module customised for the Slurm job manager

    screen
    module load bpipe/0.9.9.3.slurm

Start running GLOE-Pipe (NB: **the indirect mode is the default**)

    bpipe run gloeseq.pipeline.groovy rawdata/*.fastq.gz

## Notes
- To use the direct mode, the main pipeline needs to be changed in the `gloeseq.pipeline.groovy` file
- The `REs` folder contains bed files that include information on the break sites produced by the restriction endonucleases used in the GLOE-seq article
- For the yeast genome, the break sites that overlap with transposons, telomeres and mtDNA are excluded
